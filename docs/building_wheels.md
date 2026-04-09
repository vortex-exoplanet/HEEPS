# Building HEEPS and PyPROPER3 Wheels

## Prerequisites

- Python 3.10+ installed
- `poetry` installed (`pip install poetry`)
- `wheel` and `setuptools` installed (`pip install wheel setuptools`)

## Building the HEEPS wheel

From the repo root (`D:\Repos\HEEPS`):

```bash
python -m poetry build
```

This produces two files in `dist/`:
- `heeps-<version>-py3-none-any.whl` (wheel)
- `heeps-<version>.tar.gz` (sdist)

The version is read from `heeps/__init__.py` (`__version__`) and must match the version in `pyproject.toml`.

### What the HEEPS wheel contains

All Python modules under `heeps/`, including the `util/` subpackage. No data files (FITS, config, etc.) are bundled — those are downloaded at runtime via Google Drive.

## Building the PyPROPER3 3.3.4 wheel

PROPER 3.3.4 is **not available on PyPI** (only 3.2.1 is there), so we build it from the SourceForge zip.

```bash
# Extract the zip
mkdir /tmp/proper_build && cd /tmp/proper_build
unzip /path/to/proper_v3.3.4_python.zip
cd proper_v3.3.4_python

# Build the wheel
python setup.py bdist_wheel
```

The wheel appears in `dist/`:
- `PyPROPER3-3.3.4-py3-none-any.whl`

### Gotcha: `wheel` must be installed

If you get `error: invalid command 'bdist_wheel'`, install the `wheel` package first:

```bash
pip install wheel
```

### Gotcha: the wheel is much smaller than the zip (~160KB vs ~7MB)

This is expected. The zip contains a 7MB PDF manual (`PROPER_manual_v3.3.4.pdf`) at the top level, outside the `proper/` package directory. The `setup.py` `package_data` only bundles files inside `proper/`, so the PDF is excluded from the wheel. The `.fits` data file (72KB, inside `proper/`) is included.

## Hosting on a private PyPI server

Both wheels need to be uploaded to your private PyPI server:

1. `dist/heeps-<version>-py3-none-any.whl`
2. `PyPROPER3-3.3.4-py3-none-any.whl`

HEEPS declares `PyPROPER3>=3.3.3` as a dependency sourced from the university wheel server. The `pyproject.toml` includes a `[[tool.poetry.source]]` entry pointing to `https://scopesim.univie.ac.at/wheels/`.

For Poetry users, this is automatic. For pip users, pass the extra index:

```bash
pip install heeps --extra-index-url https://scopesim.univie.ac.at/wheels/
```

Using `--extra-index-url` (not `--index-url`) lets pip fall back to public PyPI for other dependencies like `vip-hci`, `numpy`, etc.

## Testing locally (without a PyPI server)

```bash
python -m uv venv .venv
.venv\Scripts\activate                # Windows
# source .venv/bin/activate           # Linux/macOS

python -m uv pip install /path/to/PyPROPER3-3.3.4-py3-none-any.whl
python -m uv pip install -e /path/to/HEEPS    # editable install for dev
# or: python -m uv pip install /path/to/HEEPS/dist/heeps-1.0.0-py3-none-any.whl
```

## Gotchas summary

| Issue | Cause | Fix |
|---|---|---|
| `bdist_wheel` command not found | `wheel` package not installed | `pip install wheel` |
| PROPER wheel is tiny vs zip | PDF manual excluded from wheel | Expected — not a problem |
| `pip install heeps` can't find PyPROPER3>=3.3.3 | Only 3.2.1 on public PyPI | Host 3.3.4 wheel on private server, use `--extra-index-url` |
| `uv` command not found on Windows | Installed to user Scripts dir, not on PATH | Use `python -m uv` instead |
| `heeps.util` import fails | `heeps/util/__init__.py` was missing | Already fixed — file was added |
| Version mismatch | `__init__.py` and `pyproject.toml` versions can drift | Keep both in sync when bumping |
