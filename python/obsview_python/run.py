#!/usr/bin/env python3
"""Thin executable wrapper.

Usage (from the project root, one level above the obsview/ package):

This exists because obsview/ uses relative imports internally (e.g.
`from . import config`), so its main() must be invoked either via
`python -m obsview.main` or through a wrapper like this one - not by
running obsview/main.py directly.
"""
from obsview.main import main

if __name__ == "__main__":
    main()
