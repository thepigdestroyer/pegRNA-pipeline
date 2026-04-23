"""Thin, testable wrapper around :mod:`paramiko`'s SSH client.

Design
------
The real paramiko client is obtained via a ``client_factory`` callable so
tests can inject a :class:`unittest.mock.MagicMock` without any network
activity. A typical production call site is::

    ssh = SSHClient(
        hostname="login.hpc.imperial.ac.uk",
        username="mj",
        password_provider=getpass.getpass,
    )
    with ssh:
        stdout, stderr, rc = ssh.exec("qsub job.pbs")

The class exposes synchronous ``exec``, SFTP ``upload`` / ``download``,
and ``reconnect`` for resilient polling loops.
"""

from __future__ import annotations

import socket
from typing import Any, Callable, Optional, Tuple


class SSHError(RuntimeError):
    """Base class for all :mod:`ssh_client` failures."""


class SSHAuthError(SSHError):
    """Authentication was rejected by the remote host."""


class SSHConnectError(SSHError):
    """Could not reach the remote host (DNS, TCP, banner timeout, etc.)."""


def _default_client_factory():
    """Return a fresh :class:`paramiko.SSHClient` (deferred import so
    the module stays importable in environments without paramiko)."""
    import paramiko

    return paramiko.SSHClient()


class SSHClient:
    """Password-authenticated SSH client wrapper.

    Parameters
    ----------
    hostname, username:
        Remote host + login name.
    password_provider:
        Zero-arg callable returning the password. Called at ``connect``
        time so callers can re-prompt on reconnect.
    timeout:
        TCP connect timeout in seconds (default 30).
    keepalive:
        Server-side keepalive interval in seconds (default 60). Passed to
        the underlying transport after connect.
    client_factory:
        Callable returning an object that quacks like a
        :class:`paramiko.SSHClient`. Swap this in tests.
    """

    def __init__(
        self,
        *,
        hostname: str,
        username: str,
        password_provider: Callable[[], str],
        timeout: int = 30,
        keepalive: int = 60,
        client_factory: Callable[[], Any] = _default_client_factory,
    ) -> None:
        self.hostname = hostname
        self.username = username
        self.password_provider = password_provider
        self.timeout = timeout
        self.keepalive = keepalive
        self._client_factory = client_factory
        self._client: Optional[Any] = None
        self._sftp: Optional[Any] = None
        self._connected: bool = False

    # ---------------------------------------------------------------- lifecycle
    @property
    def is_connected(self) -> bool:
        return self._connected

    def connect(self) -> None:
        """Open the SSH session, applying ``AutoAddPolicy`` for the host key."""
        import paramiko  # deferred

        client = self._client_factory()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        password = self.password_provider()
        try:
            client.connect(
                hostname=self.hostname,
                username=self.username,
                password=password,
                timeout=self.timeout,
                look_for_keys=False,
                allow_agent=False,
            )
        except paramiko.AuthenticationException as exc:
            raise SSHAuthError(str(exc)) from exc
        except (socket.gaierror, socket.timeout, OSError, paramiko.SSHException) as exc:
            raise SSHConnectError(str(exc)) from exc

        # Best-effort keepalive; tests with MagicMock won't have get_transport
        # returning a real object, so silently skip if the attribute is missing.
        try:
            transport = client.get_transport()
            if transport is not None:
                transport.set_keepalive(self.keepalive)
        except Exception:  # pragma: no cover - defensive
            pass

        self._client = client
        self._connected = True

    def disconnect(self) -> None:
        """Close the SFTP channel (if open) and the SSH transport."""
        if self._sftp is not None:
            try:
                self._sftp.close()
            except Exception:  # pragma: no cover - defensive
                pass
            self._sftp = None
        if self._client is not None:
            try:
                self._client.close()
            except Exception:  # pragma: no cover - defensive
                pass
        self._connected = False

    def reconnect(self) -> None:
        """Force a close-and-reconnect cycle (for stale sessions)."""
        self.disconnect()
        self.connect()

    # ---------------------------------------------------------------- commands
    def exec(self, command: str) -> Tuple[str, str, int]:
        """Run ``command`` on the remote host; returns ``(stdout, stderr, rc)``.

        ``stdout`` / ``stderr`` are decoded as UTF-8 (replacing bad bytes);
        ``rc`` is the process exit code.
        """
        if not self._connected or self._client is None:
            raise SSHError("not connected; call connect() first")
        _stdin, stdout, stderr = self._client.exec_command(command)
        out = stdout.read().decode("utf-8", errors="replace")
        err = stderr.read().decode("utf-8", errors="replace")
        rc = int(stdout.channel.recv_exit_status())
        return out, err, rc

    # ---------------------------------------------------------------- SFTP
    @property
    def sftp(self):
        """Lazy SFTP subclient (opens on first access, reused thereafter)."""
        if not self._connected or self._client is None:
            raise SSHError("not connected; call connect() first")
        if self._sftp is None:
            self._sftp = self._client.open_sftp()
        return self._sftp

    def upload(self, local_path: str, remote_path: str) -> None:
        """Copy a local file to a remote path (overwrites)."""
        self.sftp.put(local_path, remote_path)

    def download(self, remote_path: str, local_path: str) -> None:
        """Copy a remote file to a local path (overwrites)."""
        self.sftp.get(remote_path, local_path)

    # ---------------------------------------------------------------- context
    def __enter__(self) -> "SSHClient":
        self.connect()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.disconnect()
