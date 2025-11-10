import os
from pathlib import Path


def _remove_empty_dirs(root: Path):
    for dirpath, dirnames, filenames in os.walk(root, topdown=False):
        p = Path(dirpath)
        if not any(p.iterdir()):
            p.rmdir()


def pytest_sessionfinish(session, exitstatus):
    test_output_dir = Path(__file__).parent.parent
    if test_output_dir.exists():
        _remove_empty_dirs(test_output_dir)
