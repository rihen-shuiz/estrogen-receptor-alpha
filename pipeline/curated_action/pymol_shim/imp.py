"""Minimal Python 3.12 compatibility shim for PyMOL 2.5's `from imp import find_module`."""
import importlib.util


def find_module(name, path=None):
    spec = importlib.util.find_spec(name, path)
    if spec is None:
        raise ImportError(name)
    return (None, spec.origin, ("", "", 0))
