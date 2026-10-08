from __future__ import annotations

import ast
import io
import sys
import unittest

from process_output import diagnostic_bytes, run_checked


class ProcessOutputTest(unittest.TestCase):
    def test_round_trip_and_cp932_output(self):
        for data in (
            b"ASCII\n",
            "path/検証/測定.txt\n".encode("utf-8"),
            "path/検証/測定.txt\n".encode("cp932"),
            b"invalid UTF-8: \xff\xfe\x80\n",
            "genuine replacement: \ufffd\n".encode("utf-8"),
            bytes(range(256)),
        ):
            with self.subTest(data=data):
                displayed = diagnostic_bytes(data)
                self.assertTrue(displayed.isascii())
                self.assertEqual(ast.literal_eval(displayed), data)
                raw = io.BytesIO()
                stream = io.TextIOWrapper(raw, encoding="cp932", errors="strict")
                print(displayed, file=stream)
                stream.flush()
                self.assertIn(b"b", raw.getvalue())

    def test_success_keeps_original_bytes(self):
        out = "検証".encode("cp932") + b"\xff"
        err = "測定".encode("utf-8")
        p = run_checked([sys.executable, "-c",
                         f"import os; os.write(1,{out!r}); os.write(2,{err!r})"])
        self.assertEqual(p.returncode, 0)
        self.assertEqual(p.stdout, out)
        self.assertEqual(p.stderr, err)

    def test_expected_failure_is_not_a_skip(self):
        p = run_checked([sys.executable, "-c", "import sys; sys.exit(2)"], 2)
        self.assertEqual(p.returncode, 2)

    def test_unexpected_exit_is_reported_without_encoding_failure(self):
        data = "検証".encode("cp932") + b"\xff"
        with self.assertRaises(AssertionError) as caught:
            run_checked([sys.executable, "-c",
                         f"import os,sys; os.write(2,{data!r}); sys.exit(7)"])
        message = str(caught.exception)
        self.assertIn("expected 0, got 7", message)
        self.assertIn(repr(data), message)
        message.encode("cp932", errors="strict")


if __name__ == "__main__":
    unittest.main()
