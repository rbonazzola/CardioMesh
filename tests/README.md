## Unit tests

Install the package with its dev extras and run pytest from the repo root:

```
pip install -e ".[dev]"
pytest tests/
```

`test_cardiac3dmesh.py` covers the `Cardiac3DMesh` loaders (VTK, pickle, and
NumPy/CSV/TXT), subpart extraction, and mesh topology (edges, adjacency,
neighbors), using the sample meshes under `sample_data/`.
