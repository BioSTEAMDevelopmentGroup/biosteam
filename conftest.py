# -*- coding: utf-8 -*-
"""
Configuration for pytest to run with and without numba JIT compiled functions.
"""
import pytest
import os

def pytest_ignore_collect(path):
    path = str(path)
    if 'setup' in path:
        return True

def pytest_addoption(parser):
    parser.addoption(
        "--disable-numba", action="store", default="1", help="my option: 0 or 1"
    )

def pytest_configure(config):
    os.environ["NUMBA_DISABLE_JIT"] = config.getoption("--disable-numba")
    os.environ["DISABLE_PREFERENCES"] = "1"
    os.environ["FILTER_WARNINGS"] = "1"
