"""Local SQLite run history at ``~/.pegrna/runs.db``.

Tracks every ``pegrna run`` invocation so ``pegrna status/list/cancel/fetch``
can operate against past runs without re-contacting HPC.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from sqlalchemy import Column, DateTime, String, create_engine
from sqlalchemy.orm import DeclarativeBase, Session


class Base(DeclarativeBase):
    pass


class Run(Base):
    __tablename__ = "runs"
    run_id = Column(String, primary_key=True)
    gene = Column(String, nullable=False, default="")
    mode = Column(String, nullable=False, default="")
    construct = Column(String, nullable=False, default="")
    status = Column(String, default="submitted")
    driver_job_id = Column(String, default="")
    remote_run_dir = Column(String, default="")
    local_reports_dir = Column(String, default="")
    created_at = Column(DateTime, default=datetime.utcnow)
    completed_at = Column(DateTime)


class RunsDB:
    """Thin wrapper over SQLAlchemy sessions for CLI use."""

    @classmethod
    def local_default(cls) -> "RunsDB":
        path = Path("~/.pegrna/runs.db").expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        return cls(path)

    def __init__(self, db_path: Path) -> None:
        self.engine = create_engine(f"sqlite:///{db_path}")
        Base.metadata.create_all(self.engine)

    def record_run(
        self,
        *,
        run_id: str,
        gene: str,
        mode: str,
        construct: str,
        driver_job_id: str = "",
        remote_run_dir: str = "",
        local_reports_dir: str = "",
        status: str = "submitted",
    ) -> Run:
        with Session(self.engine) as session:
            run = Run(
                run_id=run_id,
                gene=gene,
                mode=mode,
                construct=construct,
                driver_job_id=driver_job_id,
                remote_run_dir=remote_run_dir,
                local_reports_dir=local_reports_dir,
                status=status,
            )
            session.merge(run)
            session.commit()
            return run

    def update_status(self, run_id: str, status: str) -> None:
        with Session(self.engine) as session:
            run = session.get(Run, run_id)
            if run is None:
                raise KeyError(run_id)
            run.status = status
            if status in ("success", "failed", "cancelled"):
                run.completed_at = datetime.utcnow()
            session.commit()

    def mark_cancelled(self, run_id: str) -> None:
        self.update_status(run_id, "cancelled")

    def get(self, run_id: str) -> Run:
        with Session(self.engine) as session:
            run = session.get(Run, run_id)
            if run is None:
                raise KeyError(f"no run with id {run_id!r}")
            # Detach so caller can access after session closes
            session.expunge(run)
            return run

    def list_recent(self, limit: int = 20) -> list[Run]:
        with Session(self.engine) as session:
            runs = list(session.query(Run).order_by(Run.created_at.desc()).limit(limit))
            for r in runs:
                session.expunge(r)
            return runs

    def list_all(self) -> list[Run]:
        return self.list_recent(limit=10_000)
