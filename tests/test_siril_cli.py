"""Tests for siril_cli_command: the Siril AppImage's siril-cli runs through the AppImage."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import stat
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pipeline import siril_cli_command


def _executable(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


class TestSirilCliCommand(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.appdir = self.root / ".mount_Siril-abc"
        self.appimage = str(self.root / "Siril-1.4.5-x86_64.AppImage")
        _executable(self.appdir / "usr" / "bin" / "siril-cli")
        _executable(self.root / "usr-local-bin" / "siril-cli")

    def env(self, path_dirs, appimage=True):
        env = {"PATH": os.pathsep.join(str(d) for d in path_dirs)}
        if appimage:
            env.update(APPIMAGE=self.appimage, APPDIR=str(self.appdir))
        return mock.patch.dict(os.environ, env, clear=True)

    def test_appimage_siril_cli_runs_through_the_appimage(self):
        with self.env([self.appdir / "usr" / "bin", self.root / "usr-local-bin"]):
            self.assertEqual(siril_cli_command("siril-cli"),
                             [self.appimage, "siril-cli"])

    def test_installed_siril_cli_runs_directly(self):
        # Inside the AppImage too: a siril-cli outside $APPDIR found first
        with self.env([self.root / "usr-local-bin"]):
            self.assertEqual(siril_cli_command("siril-cli"), ["siril-cli"])

    def test_no_appimage(self):
        with self.env([self.appdir / "usr" / "bin"], appimage=False):
            self.assertEqual(siril_cli_command("siril-cli"), ["siril-cli"])

    def test_absolute_path_kept(self):
        mac = str(self.root / "usr-local-bin" / "siril-cli")
        with self.env([], appimage=False):
            self.assertEqual(siril_cli_command(mac), [mac])


if __name__ == "__main__":
    unittest.main()
