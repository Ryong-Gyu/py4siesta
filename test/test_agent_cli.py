import contextlib
import io
import json
import unittest

from py4siesta import agent_cli


class AgentCliTests(unittest.TestCase):
    def test_parser_accepts_representative_commands(self):
        parser = agent_cli.build_parser()

        args = parser.parse_args(["kpoint-bulk", "--kpoints", "2", "4", "6"])
        self.assertEqual(args.command, "kpoint-bulk")
        self.assertEqual(args.kpoints, [2, 4, 6])

        args = parser.parse_args([
            "eos-sliding",
            "--selection",
            "20-30",
            "--mode",
            "fractional",
            "--vector",
            "0.25",
            "0.50",
        ])
        self.assertEqual(args.command, "eos-sliding")
        self.assertEqual(args.vector, [[0.25, 0.5]])

    def test_missing_band_file_returns_json_error(self):
        stdout = io.StringIO()
        stderr = io.StringIO()

        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = agent_cli.main(["band", "--bands-path", "missing.bands"])

        self.assertEqual(status, 1)
        self.assertEqual(stdout.getvalue(), "")
        payload = json.loads(stderr.getvalue())
        self.assertFalse(payload["ok"])
        self.assertEqual(payload["command"], "band")
        self.assertEqual(payload["error"]["type"], "FileNotFoundError")


if __name__ == "__main__":
    unittest.main()
