"""Compile and exercise the production split routines without the GEOS fixture."""

from pathlib import Path
import re
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT / "FV_StateMod.F90").read_text()


def routine(name):
    return re.search(
        rf"^subroutine {name}\b.*?^end subroutine {name}\b",
        SOURCE,
        re.MULTILINE | re.DOTALL | re.IGNORECASE,
    ).group()


class RestartSplitsTest(unittest.TestCase):
    def test_production_routines(self):
        state_type = re.search(
            r"^  type T_FVDYCORE_STATE\b.*?^  end type T_FVDYCORE_STATE\b",
            SOURCE,
            re.MULTILINE | re.DOTALL,
        ).group()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "split_state.inc").write_text(state_type)
            (path / "split_routines.inc").write_text(
                "\n".join(
                    routine(name)
                    for name in ("FV_InitSplits", "FV_SaveSplits", "adjust_fv3_splits")
                )
            )
            for precision in ([], ["-fdefault-real-8", "-fdefault-double-8"]):
                with self.subTest(precision=precision):
                    executable = path / "test_splits"
                    compilation = subprocess.run(
                        [
                            "gfortran", "-cpp", "-ffree-line-length-none",
                            "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
                            *precision, "-I", directory,
                            str(ROOT / "tests/restart_splits_driver.F90"),
                            "-o", str(executable),
                        ],
                        cwd=path, capture_output=True, text=True,
                    )
                    self.assertEqual(
                        compilation.returncode, 0,
                        compilation.stdout + compilation.stderr,
                    )
                    result = subprocess.run(
                        [str(executable), "0"], capture_output=True, text=True
                    )
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    for invalid in range(1, 8):
                        result = subprocess.run(
                            [str(executable), str(invalid)],
                            capture_output=True, text=True,
                        )
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn("restart", result.stderr.lower())

    def test_timestep_wiring(self):
        self.assertNotIn("adjust_fv3_splits", routine("State_To_FV"))
        self.assertEqual(routine("FV_Run").count("call adjust_fv3_splits"), 1)
        self.assertIn("call FV_SaveSplits(STATE)", routine("FV_To_State"))
        self.assertNotRegex(SOURCE, r"\bbase_[knq]_split\b")
        specs = (ROOT / "DynCore_GridCompMod.F90").read_text()
        block = specs.split("if (adjust_splits) then", 1)[1].split("endif", 1)[0]
        self.assertEqual(block.count("RESTART=MAPL_RestartOptional"), 4)
        self.assertEqual(block.count("DEFAULT=-1.0"), 4)


if __name__ == "__main__":
    unittest.main()
