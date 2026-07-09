import contextlib
import io
import json
import unittest

from py4siesta import agent_cli
from py4siesta.operations import KPointAnalysisOperation, KPointSamplingOperation


class AgentCliTests(unittest.TestCase):
    def test_kpoint_case_names_sort_lexically_by_sampling_value(self):
        operation = KPointSamplingOperation(context=None)

        case_names = [
            operation.case_name(index, [kpoint, kpoint, kpoint])
            for index, kpoint in enumerate([1, 2, 10], start=1)
        ]

        self.assertEqual(case_names, ["001+001+001", "002+002+002", "010+010+010"])
        self.assertEqual(sorted(case_names), case_names)
        self.assertEqual(KPointAnalysisOperation._k_value_from_case_name(case_names[-1]), 10)
        self.assertEqual(KPointAnalysisOperation._k_value_from_case_name("10+10+10"), 10)

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
