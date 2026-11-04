"""Tests for SirilRunner.run_script's timeout, with a stand-in siril-cli."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import stat
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

import pipeline
from pipeline import SirilRunner


def _fake_cli(path: Path, body: str) -> str:
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return str(path)


@unittest.skipUnless(os.name == "posix", "uses a shell script as siril-cli")
class TestRunScriptTimeout(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)
        self.log: list[str] = []
        self.runner = SirilRunner(log_cb=self.log.append)

    def run_with(self, body, timeout):
        cli = _fake_cli(self.dir / "siril-cli", body)
        with mock.patch.object(pipeline, "find_siril_cli", return_value=cli):
            t = time.monotonic()
            ok = self.runner.run_script(self.dir, ["cd ."], "_t.ssf", timeout=timeout)
            return ok, time.monotonic() - t

    def test_finished_script(self):
        ok, _ = self.run_with("echo 'log: done'\n", timeout=10)
        self.assertTrue(ok)
        self.assertFalse(self.runner.timed_out)
        self.assertIn("log: done", self.runner.last_output)

    def test_hung_script_and_its_children_are_stopped(self):
        # Like the AppImage: the real work runs in a child process.
        pid_file = self.dir / "child.pid"
        ok, took = self.run_with(
            f"echo 'log: Searching in SIMBAD'\nsleep 60 &\necho $! > {pid_file}\nwait\n",
            timeout=1)
        self.assertFalse(ok)
        self.assertTrue(self.runner.timed_out)
        self.assertLess(took, 10)
        self.assertIn("log: Searching in SIMBAD", self.runner.last_output)
        self.assertTrue(any("did not finish within 1 s" in m for m in self.log))
        child = int(pid_file.read_text())
        time.sleep(0.2)
        with self.assertRaises(ProcessLookupError):
            os.kill(child, 0)

    def test_failed_script(self):
        ok, _ = self.run_with("exit 1\n", timeout=None)
        self.assertFalse(ok)
        self.assertFalse(self.runner.timed_out)


if __name__ == "__main__":
    unittest.main()
