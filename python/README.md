# Obsview Python Port
A modern Python implementation of the GEOS-ESM observation visualization tool.
Installation:
  cd python
  pip install -e .

Optionally:
  cd python
  pip install -e ".[dev]"

This adds pytest and coverage tools for testing.

Then you can use it:

from obsview import odsload, plot_map
ods = odsload('data.ods')
fig, ax = plot_map(ods)

Or run the examples:

python examples/basic_loading.py
python examples/visualization_examples.py

Or run tests:

pytest tests/ -v
