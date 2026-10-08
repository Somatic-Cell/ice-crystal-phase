"""Test-only binary subprocess capture; diagnostics remain safe on CP932 pipes.

No guess about the child's locale, no U+FFFD replacement, no exit-code masking.
Byte repr is reversible and ASCII-only; it also distinguishes literal backslashes
from escaped bytes. This does not modify the executable or any dataset bytes.
"""
from __future__ import annotations

import os
import subprocess
from collections.abc import Sequence


def diagnostic_bytes(data: bytes) -> str:
    return repr(data)


def run_checked(
    command: Sequence[str | os.PathLike[str]], expected_code: int = 0,
) -> subprocess.CompletedProcess[bytes]:
    result = subprocess.run(command, capture_output=True)
    if result.returncode != expected_code:
        raise AssertionError(
            f"expected {expected_code}, got {result.returncode}\n"
            f"stdout bytes: {diagnostic_bytes(result.stdout)}\n"
            f"stderr bytes: {diagnostic_bytes(result.stderr)}"
        )
    return result
