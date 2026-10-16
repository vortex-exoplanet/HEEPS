# Wheel Update Reference

Two wheels are hosted on Nextcloud and installed via `pip install -r requirements.txt`:

| Wheel | Type | Source |
|---|---|---|
| `heeps-X.Y.Z-py3-none-any.whl` | pure Python | this repo |
| `pyproper3-X.Y.Z-cpABC-cpABC-linux_x86_64.whl` | C extensions (platform-specific) | SourceForge zip |

---

## Update the HEEPS wheel

```bash
# 1. Bump version in heeps/__init__.py and pyproject.toml
# 2. Rebuild
python -m build --wheel --outdir ~/wheels/ /home/gorban/github/HEEPS
# 3. Upload ~/wheels/heeps-X.Y.Z-py3-none-any.whl to Nextcloud
# 4. Update the heeps line in requirements.txt with the new share URL
```

## Update the PyPROPER wheel

```bash
# 1. Download new zip from https://sourceforge.net/projects/proper-library/files/
# 2. Extract and build
unzip proper_vX.Y.Z_python.zip -d /tmp/proper_build
cd /tmp/proper_build/proper_vX.Y.Z_python
python -m build --wheel --outdir ~/wheels/
# 3. Upload ~/wheels/pyproper3-X.Y.Z-*.whl to Nextcloud
# 4. Update the pyproper line in requirements.txt with the new share URL
```

> **Note:** The PyPROPER wheel is tied to a specific Python version and Linux architecture
> (e.g. `cp312-cp312-linux_x86_64`). Rebuild on each target platform if needed.
