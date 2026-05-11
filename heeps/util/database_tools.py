
import json
import glob
import os
import pandas as pd
import astropy.io.fits as fits
from pathlib import Path
from typing import Callable, Dict, Optional
from itables import show

import ipywidgets as widgets
from IPython.display import display, clear_output
import requests, zipfile, io, pathlib

def list_contrast_grid_db(
    db_path,
    priority_cols=None,
    drop_cols=None,
    extra_drop_cols=None,
    where: Optional[Dict[str, object]] = None,
    query: Optional[str] = None,
    predicate: Optional[Callable[[pd.DataFrame], pd.Series]] = None,
    verbose=False,
):
    """
    List all Sacred experiment entries in the HEEPS contrast grid database.

    Parameters
    ----------
    db_path : str or Path
        Path to the Sacred FileStorageObserver database directory
        (the folder that contains numbered run sub-directories).
    priority_cols : list, optional
        Columns to display first. Defaults to the main HEEPS grid parameters.
    drop_cols : list, optional
        Columns to drop entirely (replaces the default drop list).
    extra_drop_cols : list, optional
        Additional columns to drop on top of the default drop list.
    where : dict, optional
        Exact-match filter applied after building the DataFrame.
        E.g. ``{'band': 'L', 'mode': 'CVC'}``
    query : str, optional
        Pandas query string applied after building the DataFrame.
        E.g. ``"magnitude <= 5 and seeing == 'Q2'"``
    predicate : callable, optional
        Callable ``(df) -> bool Series`` applied after building the DataFrame.
        E.g. ``lambda df: df['ncpa.nmodes'] > 10``
    verbose : bool
        Print progress / warnings.

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by Sacred run id, one row per completed run.

    Examples
    --------
    # All L-band CVC runs
    list_contrast_grid_db(db_path, where={'band': 'L', 'mode': 'CVC'})

    # Query string (any pandas expression)
    list_contrast_grid_db(db_path, query="magnitude <= 5 and status == 'COMPLETED'")

    # Custom predicate
    list_contrast_grid_db(db_path, predicate=lambda df: df['ncpa.nmodes'] > 10)
    """
    db_path = Path(db_path)

    # --- Defaults ---
    if priority_cols is None:
        priority_cols = [
            'band', 'mode', 'magnitude', 'seeing',
            'scao_K', 'wv_rms', 'ncpa_freq',
            'duration', 'dit', 'hfov',
            'status', 'start_time', 'stop_time',
            'ncpa.nmodes', 'ncpa.frequency', 'ncpa.lag', 'ncpa.gain_I',
        ]

    _default_drop = [
        # long derived file paths
        'f_phase', 'f_wv', 'f_cbw', 'f_scao', 'f_amp', 'f_oat', 'f_pupil',
        'dir_output', 'dir_output_psf', 'dir_current', 'dir_input',
        # scalar derived values shown via priority_cols
        'fphase_mag', 'sigLF', 'sigHF',
        # rarely interesting defaults
        'f_lyot_stop', 'tag', 'select_lyot',
        'add_cl_vort', 'add_point_err', 'add_seg', 'add_apo_drift',
        'add_phase', 'add_amp', 'nstep', 'nframes_avg', 'dec',
        'cpu_count', 'nframes',
        # npupil / hfov are derivable from band
        'npupil',
        'do_f_phase', 'do_propagation', 'do_contrast_curves', 'dry_run',
    ]
    if drop_cols is None:
        drop_cols = _default_drop
    if extra_drop_cols:
        drop_cols = list(drop_cols) + list(extra_drop_cols)

    # --- Scan run directories ---
    run_dirs = sorted(
        [d for d in db_path.iterdir() if d.is_dir() and d.name.isdigit()],
        key=lambda d: int(d.name),
    )

    if not run_dirs:
        print(f'[list_contrast_grid_db] No numbered run directories found in {db_path}')
        return pd.DataFrame()

    records = []
    for run_dir in run_dirs:
        run_id = int(run_dir.name)
        config_file = run_dir / 'config.json'
        run_file = run_dir / 'run.json'

        if not config_file.exists():
            if verbose:
                print(f'[WARN] No config.json in run {run_id}, skipping.')
            continue

        try:
            with open(config_file) as f:
                config = json.load(f)
        except Exception as e:
            if verbose:
                print(f'[WARN] Could not read config.json for run {run_id}: {e}')
            continue

        # Flatten nested ingredient dicts  (e.g. {'ncpa': {'frequency':10}} → 'ncpa.frequency':10)
        flat = {}
        for k, v in config.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    flat[f'{k}.{kk}'] = vv
            else:
                flat[k] = v
        flat['id'] = run_id

        # Augment with Sacred run metadata
        if run_file.exists():
            try:
                with open(run_file) as f:
                    run_info = json.load(f)
                flat['status'] = run_info.get('status', 'UNKNOWN')
                flat['start_time'] = run_info.get('start_time', None)
                flat['stop_time'] = run_info.get('stop_time', None)
                flat['fail_trace'] = (run_info.get('fail_trace') or [''])[0] if run_info.get('fail_trace') else None
            except Exception:
                pass

        # Flag available artifacts
        for artifact_name in ('cc_raw.fits', 'cc_adi.fits', 'cc_adi_bckg.fits', 'conf_heeps.pkl'):
            flat[f'has_{artifact_name.replace(".", "_")}'] = (run_dir / artifact_name).exists()

        records.append(flat)

    if not records:
        print(f'[list_contrast_grid_db] No valid runs found in {db_path}')
        return pd.DataFrame()

    df = pd.DataFrame(records).set_index('id').sort_index(ascending=False)

    # Drop unwanted columns (silently ignore missing ones)
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])

    # Reorder: priority columns first, then the rest alphabetically
    existing_priority = [c for c in priority_cols if c in df.columns]
    other_cols = sorted([c for c in df.columns if c not in existing_priority])
    df = df[existing_priority + other_cols]

    # --- Optional filtering (applied after column reordering) ---
    if where:
        mask = pd.Series(True, index=df.index)
        for col, val in where.items():
            if col not in df.columns:
                raise KeyError(f"[list_contrast_grid_db] filter column '{col}' not found in DataFrame.")
            mask &= (df[col] == val)
        df = df[mask]
    if predicate is not None:
        pmask = predicate(df)
        if not isinstance(pmask, pd.Series):
            raise ValueError('[list_contrast_grid_db] predicate must return a pandas Series of booleans.')
        df = df[pmask]
    if query is not None:
        df = df.query(query)

    if verbose:
        print(f'[list_contrast_grid_db] Found {len(df)} run(s) matching the filters in {db_path}')

    return df




def download_and_extract_db(
    attachment_url,
    username,
    password,
    archive_fname,
    folder_db,
    timeout=300,
    chunk_size=1 << 20,
):
    """Download a ZIP file to disk and extract it into a target folder.

    Returns a tuple: (db_path, archive_path). Raises RuntimeError with a
    user-friendly message when download or extraction fails.
    """
    try:
        response = requests.get(
            attachment_url,
            auth=(username, password),
            timeout=timeout,
            stream=True,
        )
        response.raise_for_status()
    except requests.HTTPError as e:
        if response.status_code == 401:
            raise RuntimeError("Authentication failed - wrong username or password.") from e
        if response.status_code == 403:
            raise RuntimeError("Access denied - this account cannot read the attachment.") from e
        if response.status_code == 404:
            raise RuntimeError(
                "Attachment not found - check XWIKI_SPACE, XWIKI_PAGE, ATTACHMENT_NAME."
            ) from e
        raise RuntimeError(f"HTTP error: {e}") from e
    except requests.RequestException as e:
        raise RuntimeError(f"Download failed: {e}") from e

    archive_path = pathlib.Path(archive_fname)
    try:
        with archive_path.open("wb") as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
    except OSError as e:
        raise RuntimeError(f"Could not write downloaded file: {e}") from e

    extract_dir = pathlib.Path(folder_db)
    extract_dir.mkdir(exist_ok=True)
    try:
        with zipfile.ZipFile(archive_path) as zf:
            zf.extractall(extract_dir)
    except zipfile.BadZipFile as e:
        raise RuntimeError("Downloaded file is not a valid ZIP archive.") from e

    return str(extract_dir.resolve()), str(archive_path.resolve())


def on_download(_, status_widget, user_input_widget,
                pw_input_widget, button_widget,
                attachment_url, archive_fname, folder_db):
    """Widget callback that validates inputs and initializes DB_PATH.

    Parameters
    ----------
    _ : Any
        Button click event payload from ipywidgets.
    status_widget : widgets.Output
        Output widget used for user-facing status and error messages.
    user_input_widget : widgets.Text
        Username input widget.
    pw_input_widget : widgets.Password
        Password input widget.
    button_widget : widgets.Button
        Button that triggers the download.
    attachment_url : str
        URL of the archive to download.
    archive_fname : str
        Name of the temporary zip file name
    folder_db : str
        Folder where the archive is uncompressed
    """
    global DB_PATH
    with status_widget:
        clear_output()
        username = user_input_widget.value.strip()
        password = pw_input_widget.value
        if not username or not password:
            print("❌  Please enter both username and password.")
            return

        print("⏳  Downloading and extracting database...")
        try:
            DB_PATH, archive_path = download_and_extract_db(
                attachment_url=attachment_url,
                username=username,
                password=password,
                archive_fname=archive_fname,
                folder_db=folder_db,
            )
        except RuntimeError as e:
            print(f"❌  {e}")
            return

        button_widget.disabled = True
        user_input_widget.disabled = True
        pw_input_widget.disabled = True
        print("✅  Database ready.")
        print(f"    Downloaded file: '{archive_path}'")
        print(f"    DB_PATH = '{DB_PATH}'")

        # ── hand off to your plotting notebook ───────────────────────────────────────
        # DB_PATH is already in scope if both notebooks share the same kernel (%run).
        # Otherwise write it to a temp file so any notebook in the session can read it:
        pathlib.Path(".db_path.txt").write_text(DB_PATH)
        print("Path written to .db_path.txt")
        print("In your plotting notebook, read it with:")
        print("  import pathlib")
        print("  DB_PATH = pathlib.Path('.db_path.txt').read_text().strip()")