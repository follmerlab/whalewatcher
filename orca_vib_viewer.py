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
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import numpy as np


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


# ---------- ORCA parser ----------

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


# ---------- Main application ----------

class OrcaVibViewer(tk.Tk):
    N_FRAMES = 40     # frames per full oscillation cycle
    AMPLITUDE = 0.4   # max displacement in Angstrom (scaled to mode vector)

    def __init__(self, filepath=None):
        super().__init__()
        self.title("ORCA Vibrational Mode Viewer")
        self.geometry("1100x700")
        self.resizable(True, True)

        # Data
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

        self._build_ui()

        if filepath:
            self._load_file(filepath)

    # ------ UI construction ------

    def _build_ui(self):
        # Top bar
        topbar = tk.Frame(self, bd=1, relief=tk.RIDGE, pady=4)
        topbar.pack(fill=tk.X, side=tk.TOP)
        tk.Button(topbar, text="Open ORCA Output…", command=self._open_file,
                  font=("Helvetica", 11)).pack(side=tk.LEFT, padx=8)
        self.file_label = tk.Label(topbar, text="No file loaded", anchor="w",
                                   font=("Helvetica", 10), fg="#555555")
        self.file_label.pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)

        # Main pane: left list + right canvas
        main = tk.Frame(self)
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


# ---------- Entry point ----------

if __name__ == "__main__":
    filepath = sys.argv[1] if len(sys.argv) > 1 else None
    app = OrcaVibViewer(filepath)
    app.mainloop()
