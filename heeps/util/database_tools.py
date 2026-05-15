"""
database_tools.py — Utilities for exploring and plotting HEEPS contrast curve databases.

Functions
---------
list_contrast_grid_db
    Load and display the contrast grid database as a formatted table.
build_contrast_plotter
    Build an interactive ipywidgets-based contrast curve plotter for use in Jupyter.
build_download_widget
    Build an ipywidgets panel to download the HEEPS database from a remote source.


Generated with the assistance of GitHub Copilot (Claude Sonnet 4.6)
"""

import json
import pandas as pd
import astropy.io.fits as fits
from pathlib import Path
from typing import Callable, Dict, Optional

import ipywidgets as widgets
from IPython.display import display, clear_output
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import requests, zipfile, pathlib
import datetime

__all__ = ['list_contrast_grid_db', 'build_contrast_plotter', 'build_download_widget']
__authors__ = 'Gilles Orban de Xivry'


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


def build_download_widget(attachment_url, archive_fname, folder_db):
    """Build and return a VBox widget for credentials input and DB download.

    Creates a username Text, password Password, download Button and Output
    status area, wires them up to :func:`on_download`, and returns the
    assembled ``widgets.VBox`` ready to be passed to ``display()``.

    Parameters
    ----------
    attachment_url : str
        URL of the archive to download.
    archive_fname : str
        Name of the temporary zip file.
    folder_db : str
        Folder where the archive is extracted.

    Returns
    -------
    widgets.VBox
        The assembled widget ready to display.
    """
    user_input = widgets.Text(
        placeholder="XWiki username",
        description="Username:",
        layout=widgets.Layout(width="320px"),
    )
    pw_input = widgets.Password(
        placeholder="XWiki password",
        description="Password:",
        layout=widgets.Layout(width="320px"),
    )
    btn = widgets.Button(description="Download DB", button_style="primary")
    status = widgets.Output()

    btn.on_click(lambda b: on_download(b, status, user_input, pw_input, btn,
                                        attachment_url, archive_fname, folder_db))
    return widgets.VBox([user_input, pw_input, btn, status])


def _load_cc(run_dir: Path, artifact_name: str):
    """Load a 2-row contrast-curve FITS artifact -> (sep, cc) or (None, None)."""
    fpath = run_dir / artifact_name
    if not fpath.exists():
        return None, None
    try:
        data = fits.getdata(str(fpath))
        return data[0], data[1]
    except Exception as e:
        print(f'[WARN] Could not load {fpath}: {e}')
        return None, None


def _run_label(run_id, row):
    """Short descriptive label for a run to use in plot legend."""
    band = row.get('band', '?')
    mode = row.get('mode', '?')
    mag = row.get('magnitude', '?')
    see = row.get('seeing', '?')
    return f'#{run_id} {band}/{mode} mag={mag} s={see}'


def _run_prefix(run_id, row):
    """Unique filename prefix: run{id:04d}_{band}_{mode}_mag{mag}_s{seeing}."""
    band = row.get('band', 'X')
    mode = row.get('mode', 'X')
    mag = row.get('magnitude', 0)
    see = row.get('seeing', 'X')
    try:
        mag_str = f'{float(mag):+.2f}'.replace('+', 'p').replace('-', 'm')
    except (ValueError, TypeError):
        mag_str = str(mag)
    return f'run{int(run_id):04d}_{band}_{mode}_mag{mag_str}_s{see}'


def _apply_axes_style(ax, band, log_x=True, log_y=True):
    """Apply axis style consistent with sacred_heeps_contrast_grid.py."""
    ax.grid(True)
    ax.grid(which='minor', linestyle=':')
    if log_x and log_y:
        ax.loglog()
    elif log_x:
        ax.set_xscale('log')
    elif log_y:
        ax.set_yscale('log')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    if band in ('N1', 'N2'):
        ax.set_xticks([0.06, 0.1, 0.2, 0.5, 1, 1.2])
        ax.set_xlim(0.06, 1.3)
    else:
        ax.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5])
        ax.set_xlim(0.02, 0.9)
    ax.set_ylim(1e-8, 1e-3)


def build_contrast_plotter(df: pd.DataFrame, db_path, show_status=True, completed_only=True):
    """Build and display an interactive contrast-curve explorer with ZIP export."""
    db_path = Path(db_path)

    # Shared mutable state (current selection kept in sync with last Plot call)
    state = {'sel': pd.DataFrame(), 'fig': None}

    def _unique_sorted(col):
        if col in df.columns:
            return sorted(df[col].dropna().unique())
        return []

    bands = _unique_sorted('band')
    modes = _unique_sorted('mode')
    seeings = _unique_sorted('seeing')
    mags = _unique_sorted('magnitude')

    # Selection widgets
    style = {'description_width': '60px'}
    layout_sel = widgets.Layout(width='80px')
    layout_wide = widgets.Layout(width='160px')

    w_band = widgets.SelectMultiple(
        options=bands,
        value=bands[:1] if bands else [],
        description='',
        style=style,
        layout=layout_sel,
        rows=min(5, max(len(bands), 1)),
    )
    w_mode = widgets.SelectMultiple(
        options=modes,
        value=modes[:1] if modes else [],
        description='',
        style=style,
        layout=layout_sel,
        rows=min(5, max(len(modes), 1)),
    )
    w_seeing = widgets.SelectMultiple(
        options=seeings,
        value=seeings[:1] if seeings else [],
        description='',
        style=style,
        layout=layout_sel,
        rows=min(5, max(len(seeings), 1)),
    )
    # Magnitudes: four-column checkbox grid so all values are visible at once
    mag_checks = [widgets.Checkbox(value=True, description=str(m),
                                   style={'description_width': '0px'},
                                   layout=widgets.Layout(width='70px', margin='0'))
                  for m in mags]
    def _mag_values():
        return [mags[i] for i, cb in enumerate(mag_checks) if cb.value]
    # Split into four columns
    quarter = (len(mag_checks) + 3) // 4
    w_mag_all = widgets.ToggleButton(value=True, description='All',
                                     button_style='', icon='check',
                                     layout=widgets.Layout(width='56px', height='24px'))
    def _on_mag_all(change):
        for cb in mag_checks:
            cb.value = change['new']
    w_mag_all.observe(_on_mag_all, names='value')
    col_mag1 = widgets.VBox(mag_checks[:quarter], layout=widgets.Layout(margin='0'))
    col_mag2 = widgets.VBox(mag_checks[quarter:2*quarter], layout=widgets.Layout(margin='0 0 0 4px'))
    col_mag3 = widgets.VBox(mag_checks[2*quarter:3*quarter], layout=widgets.Layout(margin='0 0 0 4px'))
    col_mag4 = widgets.VBox(mag_checks[3*quarter:], layout=widgets.Layout(margin='0 0 0 4px'))
    mag_grid = widgets.HBox([col_mag1, col_mag2, col_mag3, col_mag4])

    # Curve-type switches
    cb_style = {'description_width': 'initial'}
    w_raw = widgets.Checkbox(value=True, description='Raw CC', style=cb_style)
    w_adi = widgets.Checkbox(value=True, description='ADI (no photon noise)', style=cb_style)
    w_bckg = widgets.Checkbox(value=True, description='ADI (with photon noise)', style=cb_style)

    # Axes options
    w_logx = widgets.Checkbox(value=True, description='Log x', style=cb_style)
    w_logy = widgets.Checkbox(value=True, description='Log y', style=cb_style)
    w_legend_inside = widgets.Checkbox(value=False, description='Legend inside', style=cb_style)

    # Run-id override
    w_ids = widgets.Text(
        value='',
        placeholder='e.g. 1,3,7  (blank = use filter)',
        description='Run IDs:',
        style={'description_width': '70px'},
        layout=widgets.Layout(width='340px'),
    )

    w_button = widgets.Button(
        description='Plot',
        button_style='primary',
        icon='bar-chart',
        layout=widgets.Layout(width='100px'),
    )

    # Status filter
    all_statuses = _unique_sorted('status')
    w_status = widgets.SelectMultiple(
        options=all_statuses,
        value=['COMPLETED'] if 'COMPLETED' in all_statuses else all_statuses,
        description='Status',
        style=style,
        layout=layout_wide,
        rows=min(5, max(len(all_statuses), 1)),
    )
    status_box = widgets.VBox(
        [widgets.HTML('<b>Status</b>'), w_status],
        layout=widgets.Layout(margin='0 8px', display='flex' if show_status else 'none'),
    )

    # Export widgets
    export_artifacts = [
        ('cc_raw.fits', 'cc_raw.fits'),
        ('cc_adi.fits', 'cc_adi.fits'),
        ('cc_adi_bckg.fits', 'cc_adi_bckg.fits'),
    ]

    w_export_dir = widgets.Text(
        value=str(Path.home()),
        description='Save zip to:',
        placeholder='directory for the zip file',
        style={'description_width': '80px'},
        layout=widgets.Layout(width='400px'),
    )
    w_export_btn = widgets.Button(
        description='Export ZIP',
        button_style='warning',
        icon='download',
        layout=widgets.Layout(width='130px'),
    )
    export_out = widgets.Output()

    # Save-plot widgets
    w_save_dir = widgets.Text(
        value=str(Path.home()),
        description='Save PDF to:',
        placeholder='directory for the PDF file',
        style={'description_width': '80px'},
        layout=widgets.Layout(width='400px'),
    )
    w_save_btn = widgets.Button(
        description='Save PDF',
        button_style='success',
        icon='file-pdf-o',
        layout=widgets.Layout(width='130px'),
    )
    save_out = widgets.Output()

    out = widgets.Output()

    # Layout
    col_curves = widgets.VBox([
        widgets.HTML('<b>Curves</b>'),
        w_raw,
        w_adi,
        w_bckg,
    ], layout=widgets.Layout(margin='0 0 0 12px', padding='0'))

    col_axes = widgets.VBox([
        widgets.HTML('<b>Axes</b>'),
        w_logx,
        w_logy,
        w_legend_inside,
    ], layout=widgets.Layout(margin='0 8px 0 4px', padding='0',
                             align_items='flex-start'))

    vdivider2 = widgets.HTML('<div style="border-left:1px solid #ccc;height:100%;margin:0 2px"></div>')

    vdivider = widgets.HTML('<div style="border-left:1px solid #ccc;height:100%;margin:0 6px"></div>')

    row_select = widgets.HBox([
        widgets.VBox([widgets.HTML('<b>Band</b>'), w_band], layout=widgets.Layout(margin='0 8px')),
        widgets.VBox([widgets.HTML('<b>Mode</b>'), w_mode], layout=widgets.Layout(margin='0 8px')),
        widgets.VBox([widgets.HTML('<b>Seeing</b>'), w_seeing], layout=widgets.Layout(margin='0 8px')),
        widgets.VBox([widgets.HBox([widgets.HTML('<b>Magnitudes</b>'), w_mag_all]), mag_grid], layout=widgets.Layout(margin='0 8px')),
        status_box,
        vdivider,
        col_curves,
        vdivider2,
        col_axes,
    ])

    col_plot = widgets.VBox([
        widgets.HTML('<b>Plot</b>'),
        w_ids,
        w_button,
    ], layout=widgets.Layout(margin='0 16px 0 8px'))

    col_save = widgets.VBox([
        widgets.HTML('<b>Save plot</b>'),
        w_save_dir,
        w_save_btn,
        save_out,
    ], layout=widgets.Layout(margin='0 16px 0 8px'))

    col_export = widgets.VBox([
        widgets.HTML('<b>Export current selection</b>'),
        w_export_dir,
        w_export_btn,
        export_out,
    ], layout=widgets.Layout(margin='0 8px'))

    row_controls = widgets.HBox([col_plot, col_save, col_export])

    ui = widgets.VBox([
        widgets.HTML('<h2 style="margin:4px 0">HEEPS Contrast Curve Explorer</h2>'),
        row_select,
        widgets.HTML('<hr style="margin:6px 0">'),
        row_controls,
        out,
    ])

    def _build_selection():
        ids_text = w_ids.value.strip()
        if ids_text:
            try:
                requested_ids = [int(x.strip()) for x in ids_text.split(',') if x.strip()]
            except ValueError:
                return None, '[ERROR] Run IDs must be comma-separated integers.'
            return df.loc[[i for i in requested_ids if i in df.index]], None

        sel = df.copy()
        filters = [
            ('band', w_band.value),
            ('mode', w_mode.value),
            ('seeing', w_seeing.value),
            ('magnitude', _mag_values()),
        ]
        if show_status:
            filters.append(('status', w_status.value))
        for col, values in filters:
            if values and col in sel.columns:
                sel = sel[sel[col].isin(values)]
        if completed_only and 'status' in sel.columns:
            sel = sel[sel['status'] == 'COMPLETED']

        sort_cols = [c for c in ('band', 'mode', 'magnitude', 'seeing') if c in sel.columns]
        if sort_cols:
            sel = sel.sort_values(sort_cols)
        return sel, None

    def on_plot(_):
        with out:
            out.clear_output(wait=True)

            sel, err = _build_selection()
            if err:
                print(err)
                return
            if sel.empty:
                print('No runs match the current selection.')
                return

            state['sel'] = sel
            print(f'Plotting {len(sel)} run(s): {list(sel.index)}')

            any_raw = w_raw.value
            any_adi = w_adi.value or w_bckg.value
            if not (any_raw or any_adi):
                print('Select at least one curve type to display.')
                return

            n_panels = int(any_raw) + int(any_adi)
            fig_width = 12 if w_legend_inside.value else 15
            fig, ax_arr = plt.subplots(
                n_panels,
                1,
                figsize=(fig_width, 4.5 * n_panels),
                sharex=False,
                squeeze=False,
            )

            ax_raw = ax_arr[0, 0] if any_raw else None
            ax_adi = ax_arr[int(any_raw), 0] if any_adi else None

            if ax_raw is not None:
                ax_raw.set_ylabel('Raw contrast')
                ax_raw.set_title('Raw contrast curve')
            if ax_adi is not None:
                ax_adi.set_ylabel('5-σ ADI sensitivity (contrast)')
                ax_adi.set_title('Post-processed 5-σ contrast')

            band = 'L'
            for i, (run_id, row) in enumerate(sel.iterrows()):
                run_dir = db_path / str(run_id)
                color = f'C{i % 10}'
                label = _run_label(run_id, row)
                band = row.get('band', 'L')

                if ax_raw is not None:
                    sep, cc = _load_cc(run_dir, 'cc_raw.fits')
                    if sep is not None:
                        ax_raw.plot(
                            sep,
                            cc,
                            color=color,
                            label=label,
                            linestyle='-',
                            marker='d',
                            markevery=0.12,
                            markersize=4,
                        )

                if ax_adi is not None:
                    if w_adi.value:
                        sep, cc = _load_cc(run_dir, 'cc_adi.fits')
                        if sep is not None:
                            ax_adi.plot(
                                sep,
                                cc,
                                color=color,
                                label=label,
                                linestyle='-',
                                marker='d',
                                markevery=0.12,
                                markersize=4,
                            )
                    if w_bckg.value:
                        sep, cc = _load_cc(run_dir, 'cc_adi_bckg.fits')
                        if sep is not None:
                            ax_adi.plot(
                                sep,
                                cc,
                                color=color,
                                label=label + ' (with ph. noise)',
                                linestyle=':',
                                marker='d',
                                markevery=0.12,
                                markersize=4,
                            )

            for ax in ax_arr.flat:
                _apply_axes_style(ax, band, log_x=w_logx.value, log_y=w_logy.value)
                ax.set_xlabel('Angular separation [arcsec]')
                if w_legend_inside.value:
                    ax.legend(fontsize=7, loc='upper right')
                else:
                    ax.legend(
                        fontsize=7,
                        loc='upper left',
                        bbox_to_anchor=(1.01, 1),
                        borderaxespad=0,
                        ncol=1,
                    )

            fig.tight_layout(rect=[0, 0, 0.75, 1] if not w_legend_inside.value else None)
            state['fig'] = fig
            plt.close(fig)
            display(fig)

    def on_save_plot(_):
        with save_out:
            save_out.clear_output(wait=True)

            fig = state['fig']
            if fig is None:
                print('[WARN] Nothing to save - click Plot first.')
                return

            out_dir = Path(w_save_dir.value.strip() or Path.home())
            if not out_dir.exists():
                try:
                    out_dir.mkdir(parents=True)
                except Exception as e:
                    print(f'[ERROR] Cannot create directory {out_dir}: {e}')
                    return

            ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
            pdf_path = out_dir / f'heeps_cc_plot_{ts}.pdf'

            try:
                fig.savefig(pdf_path, dpi=150, bbox_inches='tight')
                print('✓ Saved plot')
                print(f'  -> {pdf_path}')
            except Exception as e:
                print(f'[ERROR] Could not save PDF: {e}')

    def on_export(_):
        with export_out:
            export_out.clear_output(wait=True)

            sel = state['sel']
            if sel.empty:
                print('[WARN] Nothing to export - click Plot first.')
                return

            out_dir = Path(w_export_dir.value.strip() or Path.home())
            if not out_dir.exists():
                try:
                    out_dir.mkdir(parents=True)
                except Exception as e:
                    print(f'[ERROR] Cannot create directory {out_dir}: {e}')
                    return

            ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
            zip_path = out_dir / f'heeps_cc_export_{ts}.zip'

            n_files = 0
            with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
                for run_id, row in sel.iterrows():
                    run_dir = db_path / str(run_id)
                    prefix = _run_prefix(run_id, row)
                    for src_name, suffix in export_artifacts:
                        src = run_dir / src_name
                        if src.exists():
                            zf.write(src, arcname=f'{prefix}_{suffix}')
                            n_files += 1

            print(f'✓ Exported {n_files} file(s) from {len(sel)} run(s)')
            print(f'  -> {zip_path}')

    w_button.on_click(on_plot)
    w_save_btn.on_click(on_save_plot)
    w_export_btn.on_click(on_export)
    display(ui)
    on_plot(None)