"""Tests for DanceRunBase"""

import pytest
from dance import dancerunbase


class MockDanceRunDerived(dancerunbase.DanceRunBase):
    """A mock child of the DanceRunBase class."""

    def __init__(self):
        super().__init__()

    def run(self):
        self.check_run_fatal()


def test_derived_initializes():
    try:
        d = MockDanceRunDerived()
    except Exception as e:
        pytest.fail(f"Raised exception {e}")


def test_running_twice_causes_runtime_error():
    d = MockDanceRunDerived()
    d.run()
    with pytest.raises(RuntimeError):
        d.run()
