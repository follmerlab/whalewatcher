"""
ORCA Vibrational Mode Viewer
----------------------------
Opens an ORCA frequency/opt+freq output file, lists all vibrational frequencies,
and animates the selected normal mode as a looping 3D molecular motion.

Usage:
    python3 orca_vib_viewer.py [path/to/file.out]
"""

import sys
import re
import math
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import numpy as np
import gc
from collections import defaultdict
try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


# ---------- CPK colours and covalent radii (Angstrom) ----------
ELEMENT_COLORS = {
    "H": "#FFFFFF", "C": "#404040", "N": "#3050F8", "O": "#FF0D0D",
    "F": "#90E050", "P": "#FF8000", "S": "#FFFF30", "Cl": "#1FF01F",
    "Br": "#A62929", "I": "#940094", "Mn": "#9C7AC7", "Fe": "#E06633",
    "Co": "#F090A0", "Ni": "#50D050", "Cu": "#C88033", "Zn": "#7D80B0",
    "DEFAULT": "#BEA06E",
}
# Van-der-Waals / display radii (used purely for sphere size in plot)
DISPLAY_RADII = {
    "H": 0.31, "C": 0.77, "N": 0.75, "O": 0.73, "F": 0.71,
    "P": 1.06, "S": 1.02, "Cl": 0.99, "Br": 1.14, "I": 1.33,
    "Mn": 1.19, "Fe": 1.16, "Co": 1.11, "Ni": 1.10, "Cu": 1.12,
    "Zn": 1.18, "DEFAULT": 1.0,
}
# Covalent radii for bond detection
COV_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Mn": 1.19, "Fe": 1.16, "Co": 1.11, "Ni": 1.10, "Cu": 1.12,
    "Zn": 1.18, "DEFAULT": 1.0,
}
BOND_TOLERANCE = 0.4   # Angstrom tolerance added to sum of cov. radii

GROUP_COLOR_CYCLE = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
    "#ff7f00", "#a65628", "#f781bf", "#66c2a5",
    "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f",
]


# ---------- ORCA frequency/geometry parser ----------

def parse_orca_output(path):
    """Return (atoms, coords, freqs, modes).
    atoms : list of element symbols  (len N)
    coords: np.ndarray shape (N, 3)  Angstrom, final geometry
    freqs : list of floats           (len n_modes)
    modes : np.ndarray shape (n_modes, N, 3)  mass-weighted Cartesian displacements
    """
    with open(path, "r", errors="replace") as fh:
        lines = fh.readlines()

    # ---- 1. Grab ALL geometry blocks; keep the last one ----
    atoms = []
    coords = []
    i = 0
    last_geom_start = None
    while i < len(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in lines[i]:
            last_geom_start = i
        i += 1

    if last_geom_start is None:
        raise ValueError("No CARTESIAN COORDINATES block found.")

    i = last_geom_start + 2          # skip header + dashes line
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("-"):
            break
        parts = line.split()
        if len(parts) == 4:
            atoms.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        i += 1

    if not atoms:
        raise ValueError("Could not parse atom coordinates.")

    coords = np.array(coords)
    n_atoms = len(atoms)
    n_dof = 3 * n_atoms

    # ---- 2. Vibrational frequencies ----
    freqs = []
    freq_line_re = re.compile(r"^\s*\d+:\s+(-?\d+\.\d+)\s+cm\*\*-1")
    vib_section = False
    for line in lines:
        if "VIBRATIONAL FREQUENCIES" in line:
            vib_section = True
            continue
        if vib_section:
            m = freq_line_re.match(line)
            if m:
                freqs.append(float(m.group(1)))
            elif freqs and line.strip() == "":
                break           # blank line after last frequency → done
    if not freqs:
        raise ValueError("No vibrational frequencies found.")

    n_modes = len(freqs)

    # ---- 3. Normal modes matrix ----
    # Format: printed in blocks of 6 columns, rows 0..n_dof-1
    # Header line looks like:  "     6     7     8     9    10    11"
    modes_flat = np.zeros((n_modes, n_dof))   # modes_flat[mode_idx, dof_idx]

    in_nm_section = False
    current_col_indices = []
    # Header lines look like "                  6          7    ..." (only integers after strip)
    nm_header_re = re.compile(r"^(\d+)(\s+\d+)+$")

    for line in lines:
        if "NORMAL MODES" in line and "--------" not in line:
            in_nm_section = True
            continue
        if in_nm_section:
            if "IR SPECTRUM" in line or "RAMAN SPECTRUM" in line:
                break
            # Detect a column-header line (only integers, no decimals)
            stripped = line.strip()
            if stripped and nm_header_re.match(stripped):
                try:
                    current_col_indices = [int(x) for x in stripped.split()]
                    # Filter to valid mode indices
                    current_col_indices = [c for c in current_col_indices if c < n_modes]
                except ValueError:
                    pass
                continue
            # Detect a data row:  "   ROW   val val val ..."
            parts = stripped.split()
            if len(parts) >= 2 and current_col_indices:
                try:
                    row_idx = int(parts[0])
                    vals = [float(x) for x in parts[1:]]
                    for j, col in enumerate(current_col_indices):
                        if j < len(vals) and row_idx < n_dof and col < n_modes:
                            modes_flat[col, row_idx] = vals[j]
                except ValueError:
                    pass

    # Reshape: modes[mode, atom, xyz]
    modes = modes_flat.reshape(n_modes, n_atoms, 3)

    return atoms, coords, freqs, modes


def detect_bonds(atoms, coords, tol=BOND_TOLERANCE):
    bonds = []
    n = len(atoms)
    for i in range(n):
        ri = COV_RADII.get(atoms[i], COV_RADII["DEFAULT"])
        for j in range(i + 1, n):
            rj = COV_RADII.get(atoms[j], COV_RADII["DEFAULT"])
            d = np.linalg.norm(coords[i] - coords[j])
            if d < (ri + rj + tol):
                bonds.append((i, j))
    return bonds


# ---------- Loewdin population parser ----------

class PushbackIterator:
    """Iterator wrapper allowing one-line pushback."""
    def __init__(self, iterator):
        self.iterator = iterator
        self.pushback_line = None

    def __iter__(self):
        return self

    def __next__(self):
        if self.pushback_line is not None:
            line = self.pushback_line
            self.pushback_line = None
            return line
        return next(self.iterator)

    def push(self, line):
        self.pushback_line = line


def parse_orca_loewdin_populations_streaming(filename, chunk_size=100):
    """Stream-parse ORCA Loewdin MO populations from a .pop.log.

    Returns dict with 'spin_up' and/or 'spin_down' DataFrames.
    Rows = orbital labels (e.g. '0Cu_3dxy'), columns = MO_N.
    DataFrame.attrs carries mo_numbers, mo_energies, mo_occupations.
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for Loewdin parsing")
    results = {}
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        file_iter = PushbackIterator(iter(f))
        for line in file_iter:
            if 'SPIN UP' in line and line.strip().startswith('SPIN'):
                results['spin_up'] = _parse_loewdin_section(file_iter, chunk_size)
            elif 'SPIN DOWN' in line and line.strip().startswith('SPIN'):
                results['spin_down'] = _parse_loewdin_section(file_iter, chunk_size)
    return results


def _parse_loewdin_section(file_handle, chunk_size=100):
    all_chunks, all_mo_nums, all_energies, all_occs = [], [], [], []
    while True:
        chunk_data, spin_line = _parse_column_block(file_handle)
        if chunk_data is None:
            if spin_line is not None:
                file_handle.push(spin_line)
            break
        df_chunk, mo_nums, energies, occs = chunk_data
        if df_chunk is not None and not df_chunk.empty:
            all_chunks.append(df_chunk)
            all_mo_nums.extend(mo_nums)
            all_energies.extend(energies)
            all_occs.extend(occs)
        gc.collect()
        if spin_line is not None:
            file_handle.push(spin_line)
            break
        if df_chunk is None or df_chunk.empty:
            break
    if not all_chunks:
        return None
    full_df = pd.concat(all_chunks, axis=1)
    full_df.attrs['mo_numbers'] = all_mo_nums
    full_df.attrs['mo_energies'] = all_energies
    full_df.attrs['mo_occupations'] = all_occs
    return full_df


def _parse_column_block(file_handle):
    """Parse one block of MO columns from the Loewdin table.

    Returns ((DataFrame, mo_numbers, mo_energies, mo_occupations), spin_line).
    spin_line is set when parsing hits the next SPIN section header.
    """
    mo_numbers, mo_energies, mo_occupations = [], [], []
    orbital_data = {}
    header_lines = []
    in_data = False
    spin_line = None

    for line in file_handle:
        stripped = line.strip()
        if ('SPIN UP' in line or 'SPIN DOWN' in line) and stripped.startswith('SPIN'):
            spin_line = line
            break

        if not stripped:
            if in_data:
                break
            continue

        if 'LOEWDIN REDUCED ORBITAL' in line or 'THRESHOLD' in line:
            if in_data:
                break
            continue

        if not in_data:
            header_lines.append(line)

        if '--------' in line:
            in_data = True
            if len(header_lines) >= 4:
                for offset in [4, 3, 2]:
                    if len(header_lines) >= offset:
                        try:
                            vals = header_lines[-offset].split()
                            nums = [int(x) for x in vals]
                            if all(0 <= n < 10000 for n in nums):
                                mo_numbers = nums
                                break
                        except ValueError:
                            continue
                for offset in [3, 2]:
                    if len(header_lines) >= offset:
                        try:
                            vals = header_lines[-offset].split()
                            energies = [float(x) for x in vals]
                            if len(energies) == len(mo_numbers):
                                mo_energies = energies
                                break
                        except ValueError:
                            continue
                for offset in [2, 1]:
                    if len(header_lines) >= offset:
                        try:
                            vals = header_lines[-offset].split()
                            occs = [float(x) for x in vals]
                            if len(occs) == len(mo_numbers) and all(0 <= o <= 2 for o in occs):
                                mo_occupations = occs
                                break
                        except ValueError:
                            continue
            continue

        if in_data:
            parts = line.split()
            if len(parts) < 3 or '---' in line:
                continue
            try:
                label = f"{parts[0]}_{parts[1]}"
                pops = [float(x) for x in parts[2:]]
                if len(pops) == len(mo_numbers):
                    orbital_data[label] = pops
            except (ValueError, IndexError):
                continue

    if orbital_data and mo_numbers:
        df = pd.DataFrame.from_dict(
            orbital_data, orient='index',
            columns=[f"MO_{n}" for n in mo_numbers]
        )
        return (df, mo_numbers, mo_energies, mo_occupations), spin_line
    return None, spin_line


# ---------- Main application ----------

class OrcaVibViewer(tk.Tk):
    N_FRAMES = 40     # frames per full oscillation cycle
    AMPLITUDE = 0.4   # max displacement in Angstrom (scaled to mode vector)

    def __init__(self, filepath=None):
        super().__init__()
        self.title("ORCA Vibrational Mode Viewer")
        self.geometry("1300x760")
        self.resizable(True, True)

        # Vibrational modes data
        self.atoms = None
        self.coords = None
        self.freqs = None
        self.modes = None
        self.bonds = None
        self.anim = None
        self.current_mode = None
        self._triad_artists = []
        self._show_triad = True
        self._view_centre = None
        self._view_half_span = None

        # Loewdin orbital analysis data
        self.loewdin_data = None
        self._avail_orbitals = []
        self.orb_groups = {}    # name -> [orbital_label, ...]
        self.group_colors = {}  # name -> hex color string
        self._color_counter = 0  # monotonic; survives group deletions
        self.all_mos_var = tk.BooleanVar(value=False)
        self._last_rebuild_args = None  # cached for checkbox-driven refresh

        self._build_ui()

        if filepath:
            self._load_file(filepath)

    # ------ UI construction ------

    def _build_ui(self):
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        modes_frame = tk.Frame(self.notebook)
        self.notebook.add(modes_frame, text="Vibrational Modes")
        self._build_modes_tab(modes_frame)

        orb_frame = tk.Frame(self.notebook)
        self.notebook.add(orb_frame, text="Orbital Analysis")
        self._build_orbital_tab(orb_frame)

    def _build_modes_tab(self, parent):
        # Top bar
        topbar = tk.Frame(parent, bd=1, relief=tk.RIDGE, pady=4)
        topbar.pack(fill=tk.X, side=tk.TOP)
        tk.Button(topbar, text="Open ORCA Output…", command=self._open_file,
                  font=("Helvetica", 11)).pack(side=tk.LEFT, padx=8)
        self.file_label = tk.Label(topbar, text="No file loaded", anchor="w",
                                   font=("Helvetica", 10), fg="#555555")
        self.file_label.pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)

        # Main pane: left list + right canvas
        main = tk.Frame(parent)
        main.pack(fill=tk.BOTH, expand=True)

        # ---- Left panel: frequency list ----
        left = tk.Frame(main, width=320, bd=1, relief=tk.SUNKEN)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(6, 0), pady=6)
        left.pack_propagate(False)

        tk.Label(left, text="Vibrational Modes", font=("Helvetica", 11, "bold")
                 ).pack(pady=(6, 2))
        tk.Label(left, text="(double-click to animate)", font=("Helvetica", 9),
                 fg="#666666").pack()

        # Search / filter box
        filter_frame = tk.Frame(left)
        filter_frame.pack(fill=tk.X, padx=4, pady=4)
        tk.Label(filter_frame, text="Filter (cm⁻¹):").pack(side=tk.LEFT)
        self.filter_var = tk.StringVar()
        self.filter_var.trace_add("write", self._apply_filter)
        tk.Entry(filter_frame, textvariable=self.filter_var, width=10
                 ).pack(side=tk.LEFT, padx=2)

        # Listbox
        list_frame = tk.Frame(left)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=4, pady=2)
        scrollbar = tk.Scrollbar(list_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.listbox = tk.Listbox(list_frame, yscrollcommand=scrollbar.set,
                                  font=("Courier", 13), selectmode=tk.SINGLE,
                                  activestyle="dotbox")
        self.listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.listbox.yview)
        self.listbox.bind("<Double-1>", self._on_listbox_select)
        self.listbox.bind("<Return>", self._on_listbox_select)

        # Mode info label
        self.mode_info = tk.Label(left, text="", font=("Helvetica", 10),
                                  fg="#003366", wraplength=300, justify="left")
        self.mode_info.pack(pady=4, padx=4)

        # ---- Right panel: matplotlib figure ----
        right = tk.Frame(main)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=6, pady=6)

        self.fig = plt.Figure(figsize=(7, 5.5), dpi=100, facecolor="#1a1a2e")
        self.ax = self.fig.add_subplot(111, projection="3d",
                                       facecolor="#1a1a2e")
        self._style_axes()

        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Controls below canvas
        ctrl = tk.Frame(right)
        ctrl.pack(fill=tk.X)

        tk.Label(ctrl, text="Speed:").pack(side=tk.LEFT, padx=(0, 2))
        self.speed_var = tk.IntVar(value=40)
        speed_scale = tk.Scale(ctrl, from_=5, to=200, orient=tk.HORIZONTAL,
                               variable=self.speed_var, length=160,
                               label="", showvalue=True,
                               command=self._on_speed_change)
        speed_scale.pack(side=tk.LEFT)
        tk.Label(ctrl, text="ms/frame").pack(side=tk.LEFT)

        tk.Label(ctrl, text="  Amplitude:").pack(side=tk.LEFT, padx=(10, 2))
        self.amp_var = tk.DoubleVar(value=0.4)
        amp_scale = tk.Scale(ctrl, from_=0.05, to=1.5, resolution=0.05,
                             orient=tk.HORIZONTAL, variable=self.amp_var,
                             length=160, label="", showvalue=True,
                             command=self._on_amp_change)
        amp_scale.pack(side=tk.LEFT)
        tk.Label(ctrl, text="Å").pack(side=tk.LEFT)

        self.pause_btn = tk.Button(ctrl, text="⏸  Pause", command=self._toggle_pause,
                                   state=tk.DISABLED)
        self.pause_btn.pack(side=tk.LEFT, padx=12)

        self.axis_btn = tk.Button(ctrl, text="⊹  Axis: ON", command=self._toggle_triad,
                                  width=10)
        self.axis_btn.pack(side=tk.LEFT, padx=4)

        self._draw_placeholder()

    def _style_axes(self):
        ax = self.ax
        ax.set_facecolor("#1a1a2e")
        ax.set_axis_off()           # remove all panes, gridlines, ticks, labels
        ax.set_box_aspect([1, 1, 1])  # enforce equal physical scale on all three axes
        self.fig.tight_layout()

    # ------ File loading ------

    def _open_file(self):
        path = filedialog.askopenfilename(
            title="Select ORCA output file",
            filetypes=[("ORCA output", "*.out *.log"), ("All files", "*.*")]
        )
        if path:
            self._load_file(path)

    def _load_file(self, path):
        try:
            atoms, coords, freqs, modes = parse_orca_output(path)
        except Exception as e:
            messagebox.showerror("Parse error", str(e))
            return

        self.atoms = atoms
        self.coords = coords
        self.freqs = freqs
        self.modes = modes
        self.bonds = detect_bonds(atoms, coords)
        self.current_mode = None
        if self.anim is not None:
            self.anim.event_source.stop()
            self.anim = None

        short_name = path.split("/")[-1]
        self.file_label.config(
            text=f"{short_name}  |  {len(atoms)} atoms  |  {len(freqs)} modes"
        )
        self._populate_list(freqs)
        self._draw_placeholder()
        self.mode_info.config(text="Double-click a mode to animate")

    # ------ Frequency list ------

    def _populate_list(self, freqs=None, filter_str=""):
        if freqs is None:
            freqs = self.freqs
        self.listbox.delete(0, tk.END)
        self._list_indices = []   # map listbox row → mode index
        for i, f in enumerate(freqs):
            label = f"{i:>4d}:  {f:>10.2f} cm⁻¹"
            if filter_str:
                if filter_str not in f"{f:.2f}":
                    continue
            self.listbox.insert(tk.END, label)
            self._list_indices.append(i)
            # Colour near-zero (translation/rotation) modes grey
            if abs(f) < 5.0:
                self.listbox.itemconfig(tk.END, fg="#666666")
            elif f < 0:
                self.listbox.itemconfig(tk.END, fg="#ff6666")

    def _apply_filter(self, *_):
        self._populate_list(self.freqs, self.filter_var.get())

    def _on_listbox_select(self, event=None):
        sel = self.listbox.curselection()
        if not sel:
            return
        list_row = sel[0]
        mode_idx = self._list_indices[list_row]
        self._animate_mode(mode_idx)

    # ------ Animation ------

    def _animate_mode(self, mode_idx):
        if self.atoms is None:
            return

        # Stop any running animation
        if self.anim is not None:
            self.anim.event_source.stop()
            self.anim = None

        self.current_mode = mode_idx
        freq = self.freqs[mode_idx]
        disp = self.modes[mode_idx]   # shape (N, 3)

        # Scale displacement so max atomic displacement = amplitude
        max_norm = np.max(np.linalg.norm(disp, axis=1))
        if max_norm < 1e-10:
            messagebox.showinfo("Zero mode", f"Mode {mode_idx} has no displacement (translation/rotation).")
            return
        scale = self.amp_var.get() / max_norm

        self.mode_info.config(
            text=f"Mode {mode_idx}  |  {freq:.2f} cm⁻¹\n"
                 f"Max disp: {max_norm:.4f} (mass-wtd)"
        )

        coords0 = self.coords.copy()
        bonds = self.bonds

        # Pre-compute frames
        phases = np.sin(np.linspace(0, 2 * np.pi, self.N_FRAMES, endpoint=False))

        def _get_frame_coords(frame_idx):
            return coords0 + phases[frame_idx] * scale * disp

        # Initial draw
        self.ax.cla()
        self._style_axes()

        # Determine axis limits — equal ranges on all axes so the molecule is
        # undistorted. set_box_aspect([1,1,1]) ensures the plot box itself is cubic.
        all_coords = np.vstack([_get_frame_coords(f) for f in range(self.N_FRAMES)])
        margin = 1.5
        lo = all_coords.min(axis=0)
        hi = all_coords.max(axis=0)
        centre = (lo + hi) / 2
        half_span = max((hi - lo) / 2) + margin   # single value → uniform scale
        self.ax.set_xlim(centre[0] - half_span, centre[0] + half_span)
        self.ax.set_ylim(centre[1] - half_span, centre[1] + half_span)
        self.ax.set_zlim(centre[2] - half_span, centre[2] + half_span)
        self._view_centre = centre
        self._view_half_span = half_span

        # Atom sizes (roughly proportional to covalent radius)
        sizes = np.array([DISPLAY_RADII.get(a, DISPLAY_RADII["DEFAULT"]) * 120
                          for a in self.atoms])
        colors = [ELEMENT_COLORS.get(a, ELEMENT_COLORS["DEFAULT"]) for a in self.atoms]

        frame_coords = _get_frame_coords(0)

        # Draw bonds as lines
        bond_lines = []
        for (i, j) in bonds:
            lx = [frame_coords[i, 0], frame_coords[j, 0]]
            ly = [frame_coords[i, 1], frame_coords[j, 1]]
            lz = [frame_coords[i, 2], frame_coords[j, 2]]
            line, = self.ax.plot(lx, ly, lz, color="#888888", linewidth=1.0, zorder=1)
            bond_lines.append((i, j, line))

        # Draw atoms as scatter
        scat = self.ax.scatter(
            frame_coords[:, 0], frame_coords[:, 1], frame_coords[:, 2],
            s=sizes, c=colors, depthshade=True, zorder=2, edgecolors="#ffffff",
            linewidths=0.3
        )

        self.ax.text2D(0.5, 0.97, f"Mode {mode_idx}  —  {freq:.2f} cm⁻¹",
                       transform=self.ax.transAxes,
                       ha="center", va="top", color="#cccccc", fontsize=11)
        self._triad_artists = self._draw_triad(centre, half_span)
        for a in self._triad_artists:
            a.set_visible(self._show_triad)
        self.canvas.draw()

        # Animation update function
        def _update(frame):
            fc = _get_frame_coords(frame)
            scat._offsets3d = (fc[:, 0], fc[:, 1], fc[:, 2])
            for (i, j, line) in bond_lines:
                line.set_data_3d(
                    [fc[i, 0], fc[j, 0]],
                    [fc[i, 1], fc[j, 1]],
                    [fc[i, 2], fc[j, 2]]
                )
            return [scat] + [ln for (_, _, ln) in bond_lines]

        interval = self.speed_var.get()
        self.anim = FuncAnimation(
            self.fig, _update,
            frames=self.N_FRAMES,
            interval=interval,
            blit=False,
            repeat=True
        )
        self._paused = False
        self.pause_btn.config(text="⏸  Pause", state=tk.NORMAL)
        self.canvas.draw()

    def _draw_triad(self, centre, half_span):
        """Draw a small XYZ orientation triad in the back-bottom-left corner.
        Placed in data space so it rotates naturally with the molecule.
        Returns a list of all artists so they can be toggled."""
        # Position: back-bottom-left corner of the bounding box
        origin = centre - np.array([half_span * 0.72, half_span * 0.72, half_span * 0.72])
        arrow_len = half_span * 0.20

        axes_def = [
            (np.array([1, 0, 0]), "#ff5555", "X"),
            (np.array([0, 1, 0]), "#55dd55", "Y"),
            (np.array([0, 0, 1]), "#5588ff", "Z"),
        ]
        artists = []
        for direction, color, label in axes_def:
            tip = origin + direction * arrow_len
            line, = self.ax.plot(
                [origin[0], tip[0]], [origin[1], tip[1]], [origin[2], tip[2]],
                color=color, linewidth=2.5, zorder=10, solid_capstyle="round"
            )
            txt = self.ax.text(
                tip[0] + direction[0] * arrow_len * 0.18,
                tip[1] + direction[1] * arrow_len * 0.18,
                tip[2] + direction[2] * arrow_len * 0.18,
                label, color=color, fontsize=9, fontweight="bold", zorder=11
            )
            artists.extend([line, txt])
        return artists

    def _toggle_triad(self):
        self._show_triad = not self._show_triad
        for a in self._triad_artists:
            a.set_visible(self._show_triad)
        label = "⊹  Axis: ON" if self._show_triad else "⊹  Axis: OFF"
        self.axis_btn.config(text=label)
        self.canvas.draw_idle()

    def _draw_placeholder(self):
        self.ax.cla()
        self._style_axes()
        self.ax.text2D(0.5, 0.5, "Open an ORCA output file\nand select a vibrational mode",
                       transform=self.ax.transAxes,
                       ha="center", va="center", fontsize=13, color="#8899bb")
        self.canvas.draw()

    # ------ Controls ------

    def _toggle_pause(self):
        if self.anim is None:
            return
        if self._paused:
            self.anim.event_source.start()
            self._paused = False
            self.pause_btn.config(text="⏸  Pause")
        else:
            self.anim.event_source.stop()
            self._paused = True
            self.pause_btn.config(text="▶  Resume")

    def _on_speed_change(self, val):
        if self.anim is not None:
            self.anim.event_source.interval = int(val)

    def _on_amp_change(self, val):
        # Re-animate with new amplitude
        if self.current_mode is not None:
            self._animate_mode(self.current_mode)

    # ------ Orbital Analysis: tab layout ------

    def _build_orbital_tab(self, parent):
        if not HAS_PANDAS:
            tk.Label(
                parent,
                text="pandas is not installed.\n\nInstall with:  pip install pandas",
                font=("Helvetica", 13), fg="red"
            ).pack(expand=True)
            return

        # Top bar: file loading
        topbar = tk.Frame(parent, bd=1, relief=tk.RIDGE, pady=4)
        topbar.pack(fill=tk.X, side=tk.TOP)
        self._orb_topbar = topbar   # ref needed by _load_pop_file
        tk.Button(topbar, text="Open .pop.log…", command=self._open_pop_file,
                  font=("Helvetica", 11)).pack(side=tk.LEFT, padx=8)
        self.pop_file_label = tk.Label(topbar, text="No file loaded", anchor="w",
                                       font=("Helvetica", 10), fg="#555555")
        self.pop_file_label.pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)
        # Loading indicator — hidden until a file is parsing
        self._orb_loading_label = tk.Label(topbar, text="Loading file…",
                                           font=("Helvetica", 10, "italic"), fg="#555555")
        self._orb_progress = ttk.Progressbar(topbar, mode="indeterminate", length=160)
        # (neither is packed yet; _load_pop_file packs/unpacks them)

        # Three-column layout
        main = tk.Frame(parent)
        main.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        # ---- Column 1: Available Orbitals ----
        col1 = tk.Frame(main, bd=1, relief=tk.SUNKEN, width=200)
        col1.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 4))
        col1.pack_propagate(False)

        tk.Label(col1, text="Available Orbitals",
                 font=("Helvetica", 10, "bold")).pack(pady=(6, 2))

        filt_f = tk.Frame(col1)
        filt_f.pack(fill=tk.X, padx=4, pady=2)
        tk.Label(filt_f, text="Filter:").pack(side=tk.LEFT)
        self.orb_filter_var = tk.StringVar()
        self.orb_filter_var.trace_add("write", self._apply_orb_filter)
        tk.Entry(filt_f, textvariable=self.orb_filter_var,
                 width=12).pack(side=tk.LEFT, padx=2)

        avail_f = tk.Frame(col1)
        avail_f.pack(fill=tk.BOTH, expand=True, padx=4)
        avail_sb = tk.Scrollbar(avail_f)
        avail_sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.avail_lb = tk.Listbox(avail_f, yscrollcommand=avail_sb.set,
                                   font=("Courier", 10), selectmode=tk.EXTENDED,
                                   activestyle="dotbox", exportselection=False)
        self.avail_lb.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        avail_sb.config(command=self.avail_lb.yview)
        self.avail_lb.bind("<Double-1>", lambda e: self._add_to_group())

        tk.Button(col1, text="→ Add to Group", command=self._add_to_group,
                  font=("Helvetica", 9)).pack(fill=tk.X, padx=4, pady=(2, 4))

        # ---- Column 2: Group Management ----
        col2 = tk.Frame(main, bd=1, relief=tk.SUNKEN, width=230)
        col2.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 4))
        col2.pack_propagate(False)

        tk.Label(col2, text="Groups",
                 font=("Helvetica", 10, "bold")).pack(pady=(6, 2))

        new_f = tk.Frame(col2)
        new_f.pack(fill=tk.X, padx=4, pady=2)
        tk.Label(new_f, text="Name:").pack(side=tk.LEFT)
        self.new_grp_var = tk.StringVar(value="Group 1")
        tk.Entry(new_f, textvariable=self.new_grp_var,
                 width=9).pack(side=tk.LEFT, padx=2)
        tk.Button(new_f, text="+", width=2, command=self._new_group,
                  font=("Helvetica", 10, "bold")).pack(side=tk.LEFT)
        tk.Button(new_f, text="✕", width=2, command=self._delete_group,
                  fg="red").pack(side=tk.LEFT, padx=(2, 0))

        rename_f = tk.Frame(col2)
        rename_f.pack(fill=tk.X, padx=4, pady=(0, 2))
        tk.Button(rename_f, text="Rename selected group",
                  command=self._rename_group,
                  font=("Helvetica", 9)).pack(fill=tk.X)

        grps_f = tk.Frame(col2)
        grps_f.pack(fill=tk.X, padx=4, pady=(0, 2))
        grps_sb = tk.Scrollbar(grps_f)
        grps_sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.groups_lb = tk.Listbox(grps_f, yscrollcommand=grps_sb.set,
                                    font=("Courier", 10), selectmode=tk.SINGLE,
                                    activestyle="dotbox", height=5, exportselection=False)
        self.groups_lb.pack(side=tk.LEFT, fill=tk.X, expand=True)
        grps_sb.config(command=self.groups_lb.yview)
        self.groups_lb.bind("<<ListboxSelect>>", self._on_group_lb_select)
        self.groups_lb.bind("<Double-1>", lambda e: self._rename_group())

        tk.Label(col2, text="In selected group:",
                 font=("Helvetica", 9, "italic")).pack(pady=(4, 0), anchor="w", padx=4)

        ingrp_f = tk.Frame(col2)
        ingrp_f.pack(fill=tk.BOTH, expand=True, padx=4)
        ingrp_sb = tk.Scrollbar(ingrp_f)
        ingrp_sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.ingrp_lb = tk.Listbox(ingrp_f, yscrollcommand=ingrp_sb.set,
                                   font=("Courier", 10), selectmode=tk.EXTENDED,
                                   activestyle="dotbox", exportselection=False)
        self.ingrp_lb.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        ingrp_sb.config(command=self.ingrp_lb.yview)

        tk.Button(col2, text="← Remove Selected", command=self._remove_from_group,
                  font=("Helvetica", 9)).pack(fill=tk.X, padx=4, pady=4)

        # ---- Column 3: Plot ----
        col3 = tk.Frame(main, bd=1, relief=tk.SUNKEN)
        col3.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        ctrl_f = tk.Frame(col3)
        ctrl_f.pack(fill=tk.X, padx=6, pady=4)

        tk.Label(ctrl_f, text="Spin:").pack(side=tk.LEFT)
        self.orb_spin_var = tk.StringVar(value="up")
        ttk.Combobox(ctrl_f, textvariable=self.orb_spin_var,
                     values=["up", "down", "both"], state="readonly",
                     width=5).pack(side=tk.LEFT, padx=4)

        tk.Label(ctrl_f, text="n MOs each side:").pack(side=tk.LEFT, padx=(8, 0))
        self.n_mos_var = tk.IntVar(value=10)
        tk.Spinbox(ctrl_f, from_=1, to=100, textvariable=self.n_mos_var,
                   width=4).pack(side=tk.LEFT, padx=4)

        tk.Button(ctrl_f, text="Update Plot", command=self._update_orbital_plot,
                  font=("Helvetica", 10, "bold")).pack(side=tk.LEFT, padx=8)

        self.orb_status_label = tk.Label(ctrl_f, text="", font=("Helvetica", 9),
                                         fg="#555555")
        self.orb_status_label.pack(side=tk.LEFT)

        # Second row: gap energy display
        gap_f = tk.Frame(col3)
        gap_f.pack(fill=tk.X, padx=6, pady=(0, 4))
        tk.Label(gap_f, text="HOMO–LUMO gap:",
                 font=("Helvetica", 10, "bold")).pack(side=tk.LEFT)
        self.gap_label = tk.Label(gap_f, text="—", font=("Helvetica", 10),
                                  fg="#1a6fad")
        self.gap_label.pack(side=tk.LEFT, padx=6)

        # Display notebook: Plot tab and Table tab
        self._display_nb = ttk.Notebook(col3)
        self._display_nb.pack(fill=tk.BOTH, expand=True, padx=4, pady=(0, 4))

        # -- Plot tab --
        plot_tab = tk.Frame(self._display_nb)
        self._display_nb.add(plot_tab, text="Plot")

        self.orb_fig = plt.Figure(figsize=(7, 4.5), dpi=95)
        self.orb_ax = self.orb_fig.add_subplot(111)
        self.orb_ax.text(0.5, 0.5,
                         "Load a .pop.log, assign orbitals to groups,\nthen click Update",
                         transform=self.orb_ax.transAxes, ha="center", va="center",
                         fontsize=11, color="#888888")
        self.orb_ax.set_axis_off()
        self.orb_fig.tight_layout()

        self.orb_canvas = FigureCanvasTkAgg(self.orb_fig, master=plot_tab)
        self.orb_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.orb_canvas.draw()

        # -- Table tab --
        tbl_tab = tk.Frame(self._display_nb)
        self._display_nb.add(tbl_tab, text="Table")

        tbl_btn_f = tk.Frame(tbl_tab)
        tbl_btn_f.pack(fill=tk.X, padx=4, pady=(4, 2))
        tk.Button(tbl_btn_f, text="Copy as CSV", command=self._copy_table_csv,
                  font=("Helvetica", 9)).pack(side=tk.LEFT)
        ttk.Checkbutton(tbl_btn_f, text="Show all MOs",
                        variable=self.all_mos_var,
                        command=self._rebuild_table_refresh).pack(side=tk.LEFT, padx=10)

        self._table_frame = tk.Frame(tbl_tab)
        self._table_frame.pack(fill=tk.BOTH, expand=True, padx=4, pady=(0, 4))
        self._orb_tree = None   # built/rebuilt by _rebuild_table

    # ------ Orbital Analysis: file I/O ------

    def _open_pop_file(self):
        path = filedialog.askopenfilename(
            title="Select ORCA population log",
            filetypes=[("ORCA log", "*.log *.out"), ("All files", "*.*")]
        )
        if path:
            self._load_pop_file(path)

    def _load_pop_file(self, path):
        # Disable the open button so the user can't start a second parse
        for widget in self._orb_topbar.winfo_children():
            if isinstance(widget, tk.Button):
                widget.config(state=tk.DISABLED)

        self.pop_file_label.config(text="Loading…")
        self.orb_status_label.config(text="")

        # Show and start the indeterminate progress bar
        self._orb_loading_label.pack(side=tk.LEFT, padx=(8, 2), pady=2)
        self._orb_progress.pack(side=tk.LEFT, padx=(0, 8), pady=2)
        self._orb_progress.start(12)
        self.update_idletasks()

        result = {}

        def _parse():
            try:
                result['data'] = parse_orca_loewdin_populations_streaming(path)
            except Exception as e:
                result['error'] = str(e)

        def _check():
            if t.is_alive():
                self.after(100, _check)
                return

            # Thread done — stop the bar and re-enable the button
            self._orb_progress.stop()
            self._orb_progress.pack_forget()
            self._orb_loading_label.pack_forget()
            for widget in self._orb_topbar.winfo_children():
                if isinstance(widget, tk.Button):
                    widget.config(state=tk.NORMAL)

            if 'error' in result:
                messagebox.showerror("Parse error", result['error'])
                self.orb_status_label.config(text="Error loading file")
                self.pop_file_label.config(text="No file loaded")
                return

            data = result.get('data', {})
            if not data:
                messagebox.showwarning(
                    "No data",
                    "No SPIN UP/DOWN Loewdin MO population sections found.\n"
                    "Confirm the file contains 'LOEWDIN ORBITAL POPULATIONS PER MO'."
                )
                self.orb_status_label.config(text="No Loewdin sections found")
                self.pop_file_label.config(text="No file loaded")
                return

            self.loewdin_data = data
            all_labels = set()
            for spin_df in self.loewdin_data.values():
                if spin_df is not None:
                    all_labels.update(spin_df.index.tolist())

            self._avail_orbitals = sorted(
                all_labels,
                key=lambda s: [int(c) if c.isdigit() else c.lower()
                               for c in re.split(r'(\d+)', s)]
            )
            self._populate_orbital_list()

            short = path.replace("\\", "/").split("/")[-1]
            spins = [k.replace("spin_", "") for k in self.loewdin_data]
            self.pop_file_label.config(
                text=f"{short}  |  {len(self._avail_orbitals)} basis fns  |  spins: {', '.join(spins)}"
            )
            self.orb_status_label.config(
                text=f"Loaded {len(self._avail_orbitals)} orbital basis functions"
            )

        t = threading.Thread(target=_parse, daemon=True)
        t.start()
        self.after(100, _check)

    # ------ Orbital Analysis: orbital list ------

    def _populate_orbital_list(self, filter_str=""):
        self.avail_lb.delete(0, tk.END)
        for label in self._avail_orbitals:
            if not filter_str or filter_str.lower() in label.lower():
                self.avail_lb.insert(tk.END, label)

    def _apply_orb_filter(self, *_):
        self._populate_orbital_list(self.orb_filter_var.get())

    # ------ Orbital Analysis: group management ------

    def _current_group_name(self):
        sel = self.groups_lb.curselection()
        return self.groups_lb.get(sel[0]) if sel else None

    def _new_group(self):
        name = self.new_grp_var.get().strip()
        if not name:
            return
        if name in self.orb_groups:
            messagebox.showwarning("Duplicate", f"Group '{name}' already exists.")
            return
        idx = self._color_counter
        self._color_counter += 1
        self.orb_groups[name] = []
        self.group_colors[name] = GROUP_COLOR_CYCLE[idx % len(GROUP_COLOR_CYCLE)]
        self._refresh_groups_lb(select_name=name)
        m = re.match(r'^(.*?)(\d+)\s*$', name)
        if m:
            self.new_grp_var.set(f"{m.group(1)}{int(m.group(2)) + 1}")

    def _delete_group(self):
        name = self._current_group_name()
        if not name:
            return
        del self.orb_groups[name]
        del self.group_colors[name]
        self._refresh_groups_lb()
        self.ingrp_lb.delete(0, tk.END)

    def _rename_group(self):
        old_name = self._current_group_name()
        if not old_name:
            messagebox.showwarning("No group selected", "Select a group to rename.")
            return
        new_name = simpledialog.askstring(
            "Rename group",
            f"New name for '{old_name}':",
            initialvalue=old_name,
            parent=self
        )
        if not new_name or new_name.strip() == old_name:
            return
        new_name = new_name.strip()
        if new_name in self.orb_groups:
            messagebox.showwarning("Duplicate", f"A group named '{new_name}' already exists.")
            return
        # Rebuild dicts preserving insertion order
        self.orb_groups  = {(new_name if k == old_name else k): v
                            for k, v in self.orb_groups.items()}
        self.group_colors = {(new_name if k == old_name else k): v
                             for k, v in self.group_colors.items()}
        self._refresh_groups_lb(select_name=new_name)

    def _refresh_groups_lb(self, select_name=None):
        self.groups_lb.delete(0, tk.END)
        for i, name in enumerate(self.orb_groups):
            self.groups_lb.insert(tk.END, name)
            self.groups_lb.itemconfig(i, fg=self.group_colors.get(name, "#000000"))
        if select_name and select_name in self.orb_groups:
            idx = list(self.orb_groups.keys()).index(select_name)
            self.groups_lb.selection_set(idx)
            self._on_group_lb_select()

    def _on_group_lb_select(self, event=None):
        name = self._current_group_name()
        self.ingrp_lb.delete(0, tk.END)
        if name and name in self.orb_groups:
            for orb in self.orb_groups[name]:
                self.ingrp_lb.insert(tk.END, orb)

    def _add_to_group(self):
        name = self._current_group_name()
        if not name:
            messagebox.showwarning("No group", "Create and select a group first.")
            return
        selections = self.avail_lb.curselection()
        if not selections:
            return
        existing = set(self.orb_groups[name])
        added = 0
        for idx in selections:
            label = self.avail_lb.get(idx)
            if label not in existing:
                self.orb_groups[name].append(label)
                existing.add(label)
                added += 1
        self._on_group_lb_select()
        self.orb_status_label.config(text=f"Added {added} orbital(s) to '{name}'")

    def _remove_from_group(self):
        name = self._current_group_name()
        if not name:
            return
        to_remove = {self.ingrp_lb.get(i) for i in self.ingrp_lb.curselection()}
        self.orb_groups[name] = [o for o in self.orb_groups[name] if o not in to_remove]
        self._on_group_lb_select()

    # ------ Orbital Analysis: plotting ------

    def _update_orbital_plot(self):
        if self.loewdin_data is None:
            messagebox.showwarning("No data", "Load a .pop.log file first.")
            return
        if not self.orb_groups or all(len(v) == 0 for v in self.orb_groups.values()):
            messagebox.showwarning("No groups",
                                   "Create at least one group with orbitals assigned.")
            return

        spin = self.orb_spin_var.get()

        if spin == "both":
            results = {}
            for s in ("up", "down"):
                r = self._compute_frontier(s)
                if r is not None:
                    results[s] = r
            if not results:
                messagebox.showerror("No data", "Neither spin channel has data.")
                return

            self.orb_fig.clf()
            ncols = len(results)
            axes = self.orb_fig.subplots(1, ncols, sharey=True)
            if ncols == 1:
                axes = [axes]

            gap_parts = []
            table_data = None
            for ax, (s, r) in zip(axes, results.items()):
                mo_labels, mo_numbers, mo_energies, mo_occupations, \
                    homo_idx, lumo_idx, group_chars, frontier, df = r
                self._draw_bars(ax, mo_labels, group_chars, f"{s} spin")
                ha_to_ev = 27.2114
                homo_ev = mo_energies[homo_idx] * ha_to_ev
                lumo_ev = mo_energies[lumo_idx] * ha_to_ev
                gap_ev  = lumo_ev - homo_ev
                gap_parts.append(f"{s}: {gap_ev:.3f} eV (HOMO {homo_ev:.3f}, LUMO {lumo_ev:.3f})")
                if table_data is None:
                    table_data = (mo_labels, mo_numbers, mo_energies,
                                  mo_occupations, homo_idx, lumo_idx, group_chars, df)

            self.gap_label.config(text="   |   ".join(gap_parts))
            self.orb_fig.tight_layout()
            self.orb_canvas.draw()

            if table_data:
                self._rebuild_table(*table_data)
            self.orb_status_label.config(
                text=f"Showing both spins  |  {len(self.orb_groups)} groups"
            )
            return

        # Single spin path
        r = self._compute_frontier(spin)
        if r is None:
            messagebox.showerror("No data", f"No {spin} spin data in loaded file.")
            return
        mo_labels, mo_numbers, mo_energies, mo_occupations, \
            homo_idx, lumo_idx, group_chars, frontier, df = r

        self.orb_fig.clf()
        self.orb_ax = self.orb_fig.add_subplot(111)
        self._draw_bars(self.orb_ax, mo_labels, group_chars, f"{spin} spin")
        self.orb_fig.tight_layout()
        self.orb_canvas.draw()

        ha_to_ev = 27.2114
        homo_ev = mo_energies[homo_idx] * ha_to_ev
        lumo_ev = mo_energies[lumo_idx] * ha_to_ev
        gap_ev  = lumo_ev - homo_ev
        self.gap_label.config(
            text=f"{gap_ev:.3f} eV   (HOMO {homo_ev:.3f} eV,  LUMO {lumo_ev:.3f} eV)"
        )
        self.orb_status_label.config(
            text=f"{len(frontier)} frontier MOs  |  {len(group_chars)} groups"
        )
        self._rebuild_table(mo_labels, mo_numbers, mo_energies, mo_occupations,
                            homo_idx, lumo_idx, group_chars, df)

    def _compute_frontier(self, spin):
        """Return frontier MO data for one spin channel, or None if unavailable."""
        spin_key = f'spin_{spin}'
        if spin_key not in self.loewdin_data or self.loewdin_data[spin_key] is None:
            return None
        df = self.loewdin_data[spin_key]
        mo_numbers    = df.attrs['mo_numbers']
        mo_energies   = df.attrs['mo_energies']
        mo_occupations = df.attrs['mo_occupations']

        homo_idx = lumo_idx = None
        for i, occ in enumerate(mo_occupations):
            if occ > 0.5:
                homo_idx = i
            elif homo_idx is not None and lumo_idx is None:
                lumo_idx = i
                break
        if homo_idx is None:
            return None
        if lumo_idx is None:
            lumo_idx = len(mo_occupations) - 1

        n_mos  = self.n_mos_var.get()
        start  = max(0, homo_idx - n_mos + 1)
        end    = min(len(mo_numbers), lumo_idx + n_mos)
        frontier = list(range(start, end))

        mo_labels = []
        for i in frontier:
            if i < homo_idx:
                mo_labels.append(f'HOMO-{homo_idx - i}')
            elif i == homo_idx:
                mo_labels.append('HOMO')
            elif i == lumo_idx:
                mo_labels.append('LUMO')
            else:
                mo_labels.append(f'LUMO+{i - lumo_idx}')

        frontier_cols = [f"MO_{mo_numbers[i]}" for i in frontier]
        group_chars = {}
        for gname, gorbs in self.orb_groups.items():
            present = [o for o in gorbs if o in df.index]
            group_chars[gname] = (
                df.loc[present, frontier_cols].sum(axis=0).values
                if present else np.zeros(len(frontier))
            )

        return mo_labels, mo_numbers, mo_energies, mo_occupations, \
               homo_idx, lumo_idx, group_chars, frontier, df

    def _draw_bars(self, ax, mo_labels, group_chars, title):
        """Draw a stacked bar chart onto ax."""
        x      = np.arange(len(mo_labels))
        bottom = np.zeros(len(mo_labels))
        for gname, gchar in group_chars.items():
            color = self.group_colors.get(gname, "#999999")
            ax.bar(x, gchar, bottom=bottom, label=gname,
                   color=color, alpha=0.85, edgecolor="black", linewidth=0.4)
            bottom += gchar
        if 'HOMO' in mo_labels:
            homo_pos = mo_labels.index('HOMO')
            ax.axvline(x=homo_pos + 0.5, color='red', linestyle='--',
                       linewidth=1.5, alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(mo_labels, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel("Orbital Character (%)", fontsize=10)
        ax.set_title(f"Frontier Orbital Character  ({title})", fontsize=11)
        ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_xlim(-0.5, len(mo_labels) - 0.5)

    # ------ Orbital Analysis: table view ------

    def _rebuild_table(self, mo_labels, mo_numbers, mo_energies, mo_occupations,
                       homo_idx, lumo_idx, group_chars, df=None):
        """Rebuild the Treeview table. Respects self.all_mos_var for full-range mode."""
        # Cache args so the checkbox can trigger a refresh without a full re-plot
        self._last_rebuild_args = (mo_labels, mo_numbers, mo_energies, mo_occupations,
                                   homo_idx, lumo_idx, group_chars, df)

        ha_to_ev = 27.2114

        # Decide which MO indices to show
        if self.all_mos_var.get() and df is not None:
            indices = list(range(len(mo_numbers)))
            # Build labels: HOMO/LUMO for those two, MO number for everything else
            display_labels = []
            for i in indices:
                if i == homo_idx:
                    display_labels.append("HOMO")
                elif i == lumo_idx:
                    display_labels.append("LUMO")
                else:
                    display_labels.append(f"MO {mo_numbers[i]}")
            # Recompute group_chars for the full MO range
            all_mo_cols = [f"MO_{mo_numbers[i]}" for i in indices]
            row_group_chars = {}
            for gname, gorbs in self.orb_groups.items():
                present = [o for o in gorbs if o in df.index]
                row_group_chars[gname] = (
                    df.loc[present, all_mo_cols].sum(axis=0).values
                    if present else np.zeros(len(indices))
                )
        else:
            # Frontier slice — reconstruct from homo/lumo so indices align with group_chars
            n_mos  = self.n_mos_var.get()
            start  = max(0, homo_idx - n_mos + 1)
            end    = min(len(mo_numbers), lumo_idx + n_mos)
            indices = list(range(start, end))
            display_labels = mo_labels
            row_group_chars = group_chars

        group_names = list(row_group_chars.keys())

        # Tear down old tree + scrollbars
        for w in self._table_frame.winfo_children():
            w.destroy()

        fixed_cols = ("Orbital", "MO#", "Energy (eV)", "Occ")
        col_ids = fixed_cols + tuple(group_names)

        vsb = ttk.Scrollbar(self._table_frame, orient=tk.VERTICAL)
        hsb = ttk.Scrollbar(self._table_frame, orient=tk.HORIZONTAL)
        tree = ttk.Treeview(
            self._table_frame,
            columns=col_ids, show="headings",
            yscrollcommand=vsb.set, xscrollcommand=hsb.set,
            selectmode="browse"
        )
        vsb.config(command=tree.yview)
        hsb.config(command=tree.xview)

        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        hsb.pack(side=tk.BOTTOM, fill=tk.X)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Column headings and widths
        tree.heading("Orbital",      text="Orbital")
        tree.heading("MO#",          text="MO #")
        tree.heading("Energy (eV)",  text="Energy (eV)")
        tree.heading("Occ",          text="Occ")
        tree.column("Orbital",     width=80,  anchor=tk.CENTER)
        tree.column("MO#",         width=55,  anchor=tk.CENTER)
        tree.column("Energy (eV)", width=100, anchor=tk.E)
        tree.column("Occ",         width=45,  anchor=tk.CENTER)
        for gname in group_names:
            tree.heading(gname, text=f"{gname} (%)")
            tree.column(gname, width=max(80, len(gname) * 8 + 20), anchor=tk.E)

        tree.tag_configure("odd",  background="#f5f5f5")
        tree.tag_configure("even", background="#ffffff")
        tree.tag_configure("homo", background="#ffe0e0")
        tree.tag_configure("lumo", background="#e0e8ff")

        self._table_rows = [col_ids]
        for j, (i, orb_label) in enumerate(zip(indices, display_labels)):
            e_ev = mo_energies[i] * ha_to_ev
            occ  = mo_occupations[i]
            group_vals = [f"{row_group_chars[g][j]:.1f}" for g in group_names]

            row = (orb_label, str(mo_numbers[i]), f"{e_ev:.4f}", f"{occ:.2f}") + tuple(group_vals)
            self._table_rows.append(row)

            if orb_label == "HOMO":
                tag = "homo"
            elif orb_label == "LUMO":
                tag = "lumo"
            else:
                tag = "odd" if j % 2 else "even"
            tree.insert("", tk.END, values=row, tags=(tag,))

        self._orb_tree = tree

    def _rebuild_table_refresh(self):
        """Re-render the table using the last cached args (called when checkbox toggles)."""
        if self._last_rebuild_args is not None:
            self._rebuild_table(*self._last_rebuild_args)

    def _copy_table_csv(self):
        if not hasattr(self, '_table_rows') or not self._table_rows:
            messagebox.showinfo("No data", "Generate the table first.")
            return
        csv_text = "\n".join(",".join(str(c) for c in row) for row in self._table_rows)
        self.clipboard_clear()
        self.clipboard_append(csv_text)
        self.orb_status_label.config(text="Table copied to clipboard as CSV")


# ---------- Entry point ----------

if __name__ == "__main__":
    filepath = sys.argv[1] if len(sys.argv) > 1 else None
    app = OrcaVibViewer(filepath)
    app.mainloop()
