"""Tests for the paramiko wrapper in :mod:`pegrna_pipeline.hpc.ssh_client`.

All SSH interactions here use mocker.MagicMock so the suite never
contacts a real host. Real SSH behaviour is verified manually on the
Imperial PBS cluster.
"""

from __future__ import annotations

from typing import Any

import pytest

from pegrna_pipeline.hpc.ssh_client import (
    SSHAuthError,
    SSHClient,
    SSHConnectError,
    SSHError,
)


@pytest.fixture
def fake_client(mocker):
    """Build a MagicMock that stands in for a paramiko.SSHClient instance."""
    mock = mocker.MagicMock()
    # .connect doesn't need a return value normally.
    mock.connect.return_value = None
    return mock


def _factory(returned):
    """Return a callable compatible with the ``client_factory`` DI slot."""

    def _f():
        return returned

    return _f


def test_connect_calls_paramiko_with_password(fake_client):
    ssh = SSHClient(
        hostname="login.hpc",
        username="mj",
        password_provider=lambda: "s3cret",
        client_factory=_factory(fake_client),
    )
    ssh.connect()
    fake_client.set_missing_host_key_policy.assert_called_once()
    fake_client.connect.assert_called_once()
    kwargs = fake_client.connect.call_args.kwargs
    assert kwargs["hostname"] == "login.hpc"
    assert kwargs["username"] == "mj"
    assert kwargs["password"] == "s3cret"
    assert kwargs["timeout"] == 30


def test_connect_wraps_auth_failure(fake_client):
    import paramiko

    fake_client.connect.side_effect = paramiko.AuthenticationException("bad pw")
    ssh = SSHClient(
        hostname="h",
        username="u",
        password_provider=lambda: "x",
        client_factory=_factory(fake_client),
    )
    with pytest.raises(SSHAuthError):
        ssh.connect()


def test_connect_wraps_network_failure(fake_client):
    import socket

    fake_client.connect.side_effect = socket.gaierror("down")
    ssh = SSHClient(
        hostname="h",
        username="u",
        password_provider=lambda: "x",
        client_factory=_factory(fake_client),
    )
    with pytest.raises(SSHConnectError):
        ssh.connect()


def test_exec_returns_stdout_stderr_rc(mocker, fake_client):
    stdout = mocker.MagicMock()
    stderr = mocker.MagicMock()
    stdout.read.return_value = b"hello\n"
    stderr.read.return_value = b""
    stdout.channel.recv_exit_status.return_value = 0
    fake_client.exec_command.return_value = (mocker.MagicMock(), stdout, stderr)

    ssh = SSHClient(
        hostname="h",
        username="u",
        password_provider=lambda: "x",
        client_factory=_factory(fake_client),
    )
    ssh.connect()
    out, err, rc = ssh.exec("echo hello")
    assert out == "hello\n"
    assert err == ""
    assert rc == 0
    fake_client.exec_command.assert_called_once_with("echo hello")


def test_exec_before_connect_raises(fake_client):
    ssh = SSHClient(
        hostname="h",
        username="u",
        password_provider=lambda: "x",
        client_factory=_factory(fake_client),
    )
    with pytest.raises(SSHError):
        ssh.exec("echo x")


def test_sftp_property_lazy(fake_client):
    ssh = SSHClient(
        hostname="h",
        username="u",
        password_provider=lambda: "x",
        client_factory=_factory(fake_client),
    )
    ssh.connect()
    _ = ssh.sftp  # first access
    fake_client.open_sftp.assert_called_once()
    _ = ssh.sftp  # cached
    fake_client.open_sftp.assert_called_once()


def test_upload_uses_sftp_put(fake_client):
    ssh = SSHClient(
        hostname="h",
        username="u",
        password_provider=lambda: "x",
        client_factory=_factory(fake_client),
    )
    ssh.connect()
    sftp = fake_client.open_sftp.return_value
    ssh.upload("/tmp/local.txt", "/remote/remote.txt")
    sftp.put.assert_called_once_with("/tmp/local.txt", "/remote/remote.txt")


def test_download_uses_sftp_get(fake_client):
    ssh = SSHClient(
        hostname="h",
        username="u",
        password_provider=lambda: "x",
        client_factory=_factory(fake_client),
    )
    ssh.connect()
    sftp = fake_client.open_sftp.return_value
    ssh.download("/remote/r.txt", "/tmp/l.txt")
    sftp.get.assert_called_once_with("/remote/r.txt", "/tmp/l.txt")


def test_disconnect_closes_client(fake_client):
    ssh = SSHClient(
        hostname="h",
        username="u",
        password_provider=lambda: "x",
        client_factory=_factory(fake_client),
    )
    ssh.connect()
    ssh.disconnect()
    fake_client.close.assert_called_once()


def test_reconnect_closes_and_reconnects(fake_client):
    ssh = SSHClient(
        hostname="h",
        username="u",
        password_provider=lambda: "x",
        client_factory=_factory(fake_client),
    )
    ssh.connect()
    ssh.reconnect()
    assert fake_client.close.call_count >= 1
    assert fake_client.connect.call_count == 2


def test_context_manager(fake_client):
    with SSHClient(
        hostname="h",
        username="u",
        password_provider=lambda: "x",
        client_factory=_factory(fake_client),
    ) as ssh:
        assert ssh.is_connected
    fake_client.close.assert_called_once()
