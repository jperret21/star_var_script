#!/usr/bin/env python3
"""
Seestar S30 Pro — Variable Star Finder & Photometry  (UI entry point)
======================================================================
Run from Siril: Script > Run Script > seestar_varstar_siril.py

Requirements:
  - Siril 1.4+ open, session folder set as working directory
  - pipeline.py in the same directory as this script
  - Internet access for VizieR catalog queries

Troubleshoot:  Script > Run Script > test_connections.py
"""

from __future__ import annotations

import csv
import sys
import threading
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
import tkinter as tk
from typing import Optional

# Ensure pipeline.py (in the same folder) is importable when launched from Siril
sys.path.insert(0, str(Path(__file__).parent))

from pipeline import (  # noqa: E402
    VERSION,
    SEESTAR,
    FILTER_OPTIONS,
    HAS_SIRILPY,
    HAS_REQUESTS,
    SirilRunner,
    find_siril_cli,
    read_fits_header,
    list_fits,
    stars_in_frame,
    stars_in_safe_circle,
    find_siril_user_catalogue,
    update_siril_catalogue,
    query_vsx,
    query_apass,
    run_pipeline,
    resolve_frame_dir,
)

# ─────────────────────────────────────────────────────────────────────────────
# Theme — equilux-inspired dark palette matching Siril's UI
# ─────────────────────────────────────────────────────────────────────────────

BG      = "#2d2d2d"
BG2     = "#1f1f1f"
SURFACE = "#3c3c3c"
BORDER  = "#484848"
FG      = "#dedede"
FG2     = "#9a9a9a"
BLUE    = "#5294e2"
CYAN    = "#4eb3c9"
GREEN   = "#7ab648"
RED     = "#d45c6e"
YELLOW  = "#c89030"
MONO    = ("Menlo", "Courier New", "monospace")

# ─────────────────────────────────────────────────────────────────────────────
# Dark-themed button (tk.Label based)
# macOS Aqua ignores bg/fg on tk.Button entirely — use Label + bindings instead
# ─────────────────────────────────────────────────────────────────────────────

class DarkButton(tk.Label):
    def __init__(self, parent, text, command, hover_bg=None, **kw):
        self._cmd      = command
        self._enabled  = kw.pop("state", "normal") == "normal"
        kw.setdefault("bg",     SURFACE)
        kw.setdefault("fg",     FG)
        kw.setdefault("font",   ("Helvetica", 11))
        kw.setdefault("padx",   10)
        kw.setdefault("pady",   5)
        kw.setdefault("relief", "flat")
        kw.setdefault("cursor", "hand2")
        self._cur_bg   = kw["bg"]
        self._hover_bg = hover_bg or BORDER
        super().__init__(parent, text=text, **kw)
        self.bind("<Button-1>", lambda e: self._enabled and self._cmd())
        self.bind("<Enter>",    lambda e: self._enabled and tk.Label.configure(self, bg=self._hover_bg))
        self.bind("<Leave>",    lambda e: tk.Label.configure(self, bg=self._cur_bg))

    def configure(self, **kw):
        if "state" in kw:
            self._enabled = kw.pop("state") == "normal"
            tk.Label.configure(self, cursor="hand2" if self._enabled else "arrow")
        if "bg" in kw:
            self._cur_bg = kw["bg"]
        tk.Label.configure(self, **kw)

    config = configure

# ─────────────────────────────────────────────────────────────────────────────
# HMS/DMS → decimal degree helpers  (used when reading FITS OBJCTRA/OBJCTDEC)
# ─────────────────────────────────────────────────────────────────────────────

def _hms_to_deg(hms: str) -> float:
    """Convert 'HH MM SS.ss' (sexagesimal RA) to decimal degrees."""
    parts = hms.strip().split()
    if len(parts) >= 3:
        return (float(parts[0]) + float(parts[1]) / 60 + float(parts[2]) / 3600) * 15
    return float(hms)  # already decimal


def _dms_to_deg(dms: str) -> float:
    """Convert '±DD MM SS.ss' (sexagesimal Dec) to decimal degrees."""
    dms = dms.strip()
    sign = -1 if dms.startswith('-') else 1
    parts = dms.lstrip('+-').split()
    if len(parts) >= 3:
        return sign * (float(parts[0]) + float(parts[1]) / 60 + float(parts[2]) / 3600)
    return float(dms)  # already decimal


# ─────────────────────────────────────────────────────────────────────────────
# Main window
# ─────────────────────────────────────────────────────────────────────────────

class App(tk.Tk):
    def __init__(self, session_dir: Optional[Path] = None):
        super().__init__()
        self.title(f"Seestar S30 Pro — Variable Star Photometry  v{VERSION}")
        self.configure(bg=BG)
        self.minsize(980, 700)

        self.session_dir:    Optional[Path]  = session_dir
        self.field_ra:       float           = 0.0
        self.field_dec:      float           = 0.0
        self.field_wcs_hdr:  dict            = {}
        self.field_naxis1:   int             = 0
        self.field_naxis2:   int             = 0
        self.all_stars:      list[dict]      = []
        self.comp_stars:     list[dict]      = []
        self.selected_star:  Optional[dict]  = None
        self.runner:         Optional[SirilRunner] = None

        self._build_ui()
        self._apply_ttk_style()
        self.after(100, self._init_bg)

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        hdr = tk.Label(self,
                       text=f"  Seestar S30 Pro — Variable Star Finder & Photometry  v{VERSION}",
                       bg=BG, fg=BLUE,
                       font=("Helvetica", 15, "bold"), anchor="w")
        hdr.pack(fill="x", padx=12, pady=(12, 2))

        tk.Frame(self, bg=BORDER, height=1).pack(fill="x", padx=12)

        row = self._row(self)
        row.pack(fill="x", padx=12, pady=6)
        tk.Label(row, text="Session:", bg=BG, fg=FG2,
                 font=("Helvetica", 12)).pack(side="left")
        self.sess_var = tk.StringVar()
        tk.Entry(row, textvariable=self.sess_var, bg=BG2, fg=FG,
                 insertbackground=FG, font=("Helvetica", 12), width=52,
                 relief="flat", highlightthickness=1,
                 highlightbackground=BORDER,
                 highlightcolor=BLUE).pack(side="left", padx=(6, 4))
        self._button(row, "Browse…",  self._browse).pack(side="left", padx=3)
        self._button(row, "Load",     self._load_session).pack(side="left")

        info = self._row(self)
        info.pack(fill="x", padx=12, pady=(0, 4))
        self.lbl_field = tk.Label(info, text="Field: —", bg=BG, fg=YELLOW,
                                  font=("Helvetica", 11))
        self.lbl_field.pack(side="left")
        self.lbl_siril = tk.Label(info, text="", bg=BG, fg=FG2,
                                  font=("Helvetica", 10))
        self.lbl_siril.pack(side="right")

        tk.Frame(self, bg=BORDER, height=1).pack(fill="x", padx=12)

        paned = tk.PanedWindow(self, orient="horizontal",
                               bg=BORDER, sashwidth=4, sashrelief="flat")
        paned.pack(fill="both", expand=True, padx=12, pady=8)

        # Left panel — star table
        left = tk.Frame(paned, bg=BG)
        paned.add(left, minsize=500)

        tbl_hdr = self._row(left)
        tbl_hdr.pack(fill="x", pady=(0, 4))
        tk.Label(tbl_hdr, text="Variable stars in field",
                 bg=BG, fg=BLUE, font=("Helvetica", 12, "bold")).pack(side="left")
        self.lbl_count = tk.Label(tbl_hdr, text="0 stars", bg=BG, fg=FG2,
                                  font=("Helvetica", 11))
        self.lbl_count.pack(side="right")
        self._button(tbl_hdr, "↻  Query VSX catalog",
                     self._fetch_vsx).pack(side="right", padx=6)

        tree_frame = tk.Frame(left, bg=BG)
        tree_frame.pack(fill="both", expand=True)
        vsb = ttk.Scrollbar(tree_frame, orient="vertical")
        vsb.pack(side="right", fill="y")
        cols = ("name", "type", "period", "mag", "dist")
        self.tree = ttk.Treeview(tree_frame, columns=cols, show="headings",
                                 selectmode="browse", style="App.Treeview",
                                 yscrollcommand=vsb.set)
        headers = [("name","Name",230), ("type","Type",75),
                   ("period","Period",100), ("mag","Mag",60), ("dist","Dist°",65)]
        for col, label, w in headers:
            self.tree.heading(col, text=label,
                              command=lambda c=col: self._sort(c))
            anchor = "w" if col == "name" else "center"
            self.tree.column(col, width=w, anchor=anchor, stretch=(col=="name"))
        self.tree.pack(side="left", fill="both", expand=True)
        vsb.configure(command=self.tree.yview)
        self.tree.bind("<<TreeviewSelect>>", self._on_select)

        filt_row = self._row(left)
        filt_row.pack(fill="x", pady=(4, 0))
        tk.Label(filt_row, text="Filter:", bg=BG, fg=FG2,
                 font=("Helvetica", 11)).pack(side="left")
        self.filter_var = tk.StringVar()
        self.filter_var.trace_add("write", lambda *_: self._refresh_table())
        tk.Entry(filt_row, textvariable=self.filter_var, bg=BG2, fg=FG,
                 insertbackground=FG, font=("Helvetica", 11), width=18,
                 relief="flat", highlightthickness=1,
                 highlightbackground=BORDER,
                 highlightcolor=BLUE).pack(side="left", padx=6)
        tk.Label(filt_row, text="Max mag:", bg=BG, fg=FG2,
                 font=("Helvetica", 11)).pack(side="left")
        self.mag_max_var = tk.StringVar(value="14")
        self.mag_max_var.trace_add("write", lambda *_: self._refresh_table())
        tk.Entry(filt_row, textvariable=self.mag_max_var, bg=BG2, fg=FG,
                 insertbackground=FG, font=("Helvetica", 11), width=5,
                 relief="flat", highlightthickness=1,
                 highlightbackground=BORDER,
                 highlightcolor=BLUE).pack(side="left", padx=4)

        # Right panel
        right = tk.Frame(paned, bg=BG)
        paned.add(right, minsize=380)

        card = tk.LabelFrame(right, text=" Selected target ",
                             bg=BG, fg=BLUE, font=("Helvetica", 11, "bold"),
                             relief="solid", bd=1, padx=8, pady=6)
        card.pack(fill="x", pady=(0, 8))

        self.lbl_name = tk.Label(card, text="—", bg=BG, fg=RED,
                                 font=("Helvetica", 15, "bold"), anchor="w")
        self.lbl_name.pack(fill="x")
        self.lbl_type = tk.Label(card, text="", bg=BG, fg=CYAN,
                                 font=("Helvetica", 11), anchor="w")
        self.lbl_type.pack(fill="x")
        self.lbl_coord = tk.Label(card, text="", bg=BG, fg=FG2,
                                  font=("Helvetica", 10), anchor="w")
        self.lbl_coord.pack(fill="x")

        obs_card = tk.LabelFrame(right, text=" Observation settings ",
                                 bg=BG, fg=BLUE, font=("Helvetica", 11, "bold"),
                                 relief="solid", bd=1, padx=8, pady=6)
        obs_card.pack(fill="x", pady=(0, 6))

        flt_row = self._row(obs_card)
        flt_row.pack(fill="x", pady=(0, 4))
        tk.Label(flt_row, text="Filter:", bg=BG, fg=FG2,
                 font=("Helvetica", 11), width=9, anchor="w").pack(side="left")
        self.obs_filter_var = tk.StringVar(value="LP anti-pollution filter (CV)")
        ttk.Combobox(flt_row, textvariable=self.obs_filter_var,
                     values=list(FILTER_OPTIONS.keys()),
                     state="readonly", width=30,
                     style="Dark.TCombobox").pack(side="left", padx=(0, 4))

        step_row = self._row(obs_card)
        step_row.pack(fill="x")
        tk.Label(step_row, text="Start from:", bg=BG, fg=FG2,
                 font=("Helvetica", 11), width=9, anchor="w").pack(side="left")
        self.start_from_var = tk.StringVar(value="Full pipeline (steps 1–5)")
        ttk.Combobox(step_row, textvariable=self.start_from_var,
                     values=["Full pipeline (steps 1–5)",
                             "Apply reg + plate solve (steps 3–5)",
                             "Plate solve only (step 4–5)",
                             "Photometry only (step 5)"],
                     state="readonly", width=34,
                     style="Dark.TCombobox").pack(side="left", padx=(0, 4))

        nstars_row = self._row(obs_card)
        nstars_row.pack(fill="x", pady=(4, 0))
        tk.Label(nstars_row, text="Comp stars:", bg=BG, fg=FG2,
                 font=("Helvetica", 11), width=9, anchor="w").pack(side="left")
        self.nstars_var = tk.IntVar(value=10)
        tk.Spinbox(nstars_row, from_=3, to=19, textvariable=self.nstars_var,
                   width=5, bg=BG2, fg=FG, insertbackground=FG,
                   buttonbackground=SURFACE, font=("Helvetica", 11),
                   relief="flat", highlightthickness=1,
                   highlightbackground=BORDER,
                   highlightcolor=BLUE).pack(side="left", padx=(0, 8))
        tk.Label(nstars_row, text="(3–19, Siril max=19)", bg=BG, fg=FG2,
                 font=("Helvetica", 10)).pack(side="left")

        phot_card = tk.LabelFrame(right, text=" Photometry (setphot) ",
                                  bg=BG, fg=BLUE, font=("Helvetica", 11, "bold"),
                                  relief="solid", bd=1, padx=8, pady=6)
        phot_card.pack(fill="x", pady=(0, 6))

        self.ap_var      = tk.DoubleVar(value=10.0)
        self.inner_var   = tk.DoubleVar(value=20.0)
        self.outer_var   = tk.DoubleVar(value=30.0)
        self.dyn_var     = tk.DoubleVar(value=4.0)
        self.maxval_var  = tk.DoubleVar(value=60000.0)

        def _phot_spin(parent, var, frm, to, inc, width=7):
            return tk.Spinbox(parent, from_=frm, to=to, increment=inc,
                              textvariable=var, width=width, bg=BG2, fg=FG,
                              insertbackground=FG, buttonbackground=SURFACE,
                              font=("Helvetica", 11), relief="flat",
                              highlightthickness=1, highlightbackground=BORDER,
                              highlightcolor=BLUE)

        pr1 = self._row(phot_card)
        pr1.pack(fill="x", pady=(0, 4))
        tk.Label(pr1, text="Inner:", bg=BG, fg=FG2, font=("Helvetica", 11),
                 width=9, anchor="w").pack(side="left")
        _phot_spin(pr1, self.inner_var, 3, 300, 5).pack(side="left", padx=(0, 8))
        tk.Label(pr1, text="Outer:", bg=BG, fg=FG2, font=("Helvetica", 11),
                 anchor="w").pack(side="left")
        _phot_spin(pr1, self.outer_var, 5, 400, 5).pack(side="left", padx=(4, 0))

        pr2 = self._row(phot_card)
        pr2.pack(fill="x", pady=(0, 4))
        tk.Label(pr2, text="Dyn ratio:", bg=BG, fg=FG2, font=("Helvetica", 11),
                 width=9, anchor="w").pack(side="left")
        _phot_spin(pr2, self.dyn_var, 1.0, 5.0, 0.5).pack(side="left", padx=(0, 8))
        tk.Label(pr2, text="Aper:", bg=BG, fg=FG2, font=("Helvetica", 11),
                 anchor="w").pack(side="left")
        _phot_spin(pr2, self.ap_var, 2, 100, 1).pack(side="left", padx=(4, 0))

        pr3 = self._row(phot_card)
        pr3.pack(fill="x")
        tk.Label(pr3, text="Max val:", bg=BG, fg=FG2, font=("Helvetica", 11),
                 width=9, anchor="w").pack(side="left")
        _phot_spin(pr3, self.maxval_var, 1000, 65535, 1000, width=9).pack(
            side="left", padx=(0, 8))
        tk.Label(pr3, text="(saturation cut for comp stars)", bg=BG, fg=FG2,
                 font=("Helvetica", 10)).pack(side="left")

        cal_card = tk.LabelFrame(right, text=" Calibration frames  (optional) ",
                                 bg=BG, fg=BLUE, font=("Helvetica", 11, "bold"),
                                 relief="solid", bd=1, padx=8, pady=6)
        cal_card.pack(fill="x", pady=(0, 8))

        self.dark_var = tk.StringVar()
        self.flat_var = tk.StringVar()
        self.bias_var = tk.StringVar()
        for label, var in [("Darks:", self.dark_var),
                            ("Flats:", self.flat_var),
                            ("Bias:",  self.bias_var)]:
            row = self._row(cal_card)
            row.pack(fill="x", pady=2)
            tk.Label(row, text=label, bg=BG, fg=FG2,
                     font=("Helvetica", 11), width=7, anchor="w").pack(side="left")
            tk.Entry(row, textvariable=var, bg=BG2, fg=FG,
                     insertbackground=FG, font=("Helvetica", 10), width=22,
                     relief="flat", highlightthickness=1,
                     highlightbackground=BORDER,
                     highlightcolor=BLUE).pack(side="left", padx=(0, 4))
            self._button(row, "Browse…",
                         lambda v=var: self._browse_calib(v)).pack(side="left")

        self.run_btn = DarkButton(
            right,
            text="▶   Prepare frames & Generate light curve",
            command=self._run_pipeline,
            state="disabled",
            bg=SURFACE,
            fg=FG2,
            hover_bg="#3a7bd5",
            font=("Helvetica", 13, "bold"),
            pady=10,
        )
        self.run_btn.pack(fill="x", pady=8)

        self.prog_var = tk.IntVar()
        prog = ttk.Progressbar(right, variable=self.prog_var, maximum=100,
                               style="Green.Horizontal.TProgressbar")
        prog.pack(fill="x", pady=(0, 3))
        self.lbl_prog = tk.Label(right, text="Ready", bg=BG, fg=GREEN,
                                 font=("Helvetica", 11), anchor="w")
        self.lbl_prog.pack(fill="x")

        log_card = tk.LabelFrame(right, text=" Siril log  (select & copy with Cmd+C) ",
                                 bg=BG, fg=BLUE, font=("Helvetica", 11, "bold"),
                                 relief="solid", bd=1)
        log_card.pack(fill="both", expand=True, pady=(8, 0))

        log_frame = tk.Frame(log_card, bg=BG2)
        log_frame.pack(fill="both", expand=True, padx=4, pady=4)
        log_sb = ttk.Scrollbar(log_frame, orient="vertical")
        log_sb.pack(side="right", fill="y")
        self.log = tk.Text(log_frame, bg=BG2, fg=FG2,
                           font=(MONO[0], 10),
                           relief="flat", wrap="none",
                           state="normal",
                           height=10, yscrollcommand=log_sb.set)
        self.log.pack(side="left", fill="both", expand=True)
        log_sb.configure(command=self.log.yview)

    def _apply_ttk_style(self):
        s = ttk.Style()
        s.theme_use("clam")
        s.configure("App.Treeview",
                    background=BG2, foreground=FG,
                    fieldbackground=BG2, rowheight=23,
                    font=("Helvetica", 11))
        s.configure("App.Treeview.Heading",
                    background=SURFACE, foreground=BLUE,
                    font=("Helvetica", 11, "bold"), relief="flat")
        s.map("App.Treeview",
              background=[("selected", SURFACE)],
              foreground=[("selected", BLUE)])
        s.configure("Green.Horizontal.TProgressbar",
                    troughcolor=BG2, background=BLUE,
                    lightcolor="#3a7bd5", darkcolor="#2060c0")
        s.configure("TScrollbar",
                    background=SURFACE, troughcolor=BG2,
                    arrowcolor=FG2, bordercolor=BORDER)
        s.configure("Dark.TCombobox",
                    fieldbackground=BG2, background=SURFACE,
                    foreground=FG, arrowcolor=FG2,
                    bordercolor=BORDER, lightcolor=BG, darkcolor=BG)
        s.map("Dark.TCombobox",
              fieldbackground=[("readonly", BG2)],
              selectbackground=[("readonly", SURFACE)],
              selectforeground=[("readonly", FG)])

    def _row(self, parent):
        return tk.Frame(parent, bg=BG)

    def _button(self, parent, text, cmd):
        return DarkButton(parent, text=text, command=cmd)

    def _browse_calib(self, var: tk.StringVar):
        d = filedialog.askdirectory(title="Select folder containing calibration frames")
        if d:
            var.set(d)

    # ── Init ──────────────────────────────────────────────────────────────────

    def _init_bg(self):
        threading.Thread(target=self._connect_and_detect, daemon=True).start()

    def _connect_and_detect(self):
        self.runner = SirilRunner(log_cb=self._log)
        cli = find_siril_cli()
        if HAS_SIRILPY and self.runner._iface:
            status, color = "● Siril connected", GREEN
        elif cli:
            status, color = "● siril-cli ready", CYAN
        else:
            status, color = "✗ siril-cli not found", RED
        self.after(0, lambda: self.lbl_siril.configure(text=status, fg=color))

        if self.session_dir is None:
            # A session holds lights in either layout: lights/*.fits (legacy) or
            # Lights/<filter>/*.fits (Argos).
            def _is_session(d: Path) -> bool:
                return d.is_dir() and resolve_frame_dir(d, ["lights", "light"]) is not None

            wd = self.runner.get_working_dir()
            if wd and _is_session(wd):
                self.session_dir = wd
            else:
                for base in [Path("/Volumes/Seestar/MyWorks"),
                             Path.home() / "Argos" / "sessions",
                             Path.home() / "Desktop",
                             Path.home() / "Documents"]:
                    if not base.exists():
                        continue
                    for child in sorted(base.iterdir(), reverse=True):
                        if _is_session(child):
                            self.session_dir = child
                            break
                    if self.session_dir:
                        break

        if self.session_dir:
            self.after(0, lambda: self.sess_var.set(str(self.session_dir)))
            self.after(0, self._load_session)

    # ── Session ───────────────────────────────────────────────────────────────

    def _browse(self):
        d = filedialog.askdirectory(title="Select session folder",
                                    initialdir=str(Path.home()))
        if d:
            self.sess_var.set(d)
            self._load_session()

    def _load_session(self):
        s = self.sess_var.get().strip()
        if not s:
            return
        self.session_dir = Path(s)
        if not self.session_dir.is_dir():
            messagebox.showerror("Not found", str(self.session_dir))
            return
        threading.Thread(target=self._load_session_bg, daemon=True).start()

    def _load_session_bg(self):
        session = self.session_dir

        # Auto-fill calibration fields from the session's own calib folders
        # (Argos: Darks/<sub>/, Flats/<sub>/, Biases/<sub>/) unless already set.
        for var, names in [(self.dark_var, ["darks", "dark"]),
                           (self.flat_var, ["flats", "flat"]),
                           (self.bias_var, ["biases", "bias"])]:
            if not var.get().strip():
                d = resolve_frame_dir(session, names)
                if d:
                    self.after(0, lambda v=var, p=d: v.set(str(p)))
                    self._log(f"Calibration auto-detected: {d.name} → {d}")

        stack = None
        for cand in ["process/lights.fit", "process/result.fit"]:
            p = session / cand
            if p.exists():
                stack = p
                break
        if stack is None:
            for p in session.glob("*_og.fit"):
                stack = p
                break

        ra = dec = 0.0
        obj = ""
        n_frames = 0
        h: dict = {}
        if stack:
            h = read_fits_header(stack)
            ra       = float(h.get("RA",       h.get("CRVAL1", 0)))
            dec      = float(h.get("DEC",      h.get("CRVAL2", 0)))
            obj      = str(h.get("OBJECT",     "?"))
            n_frames = int(float(h.get("STACKCNT", 0)))
        else:
            # No stack yet — read RA/Dec from the first light frame
            lights = resolve_frame_dir(session, ["lights", "light"])
            if lights:
                first = next(iter(list_fits(lights)), None)
                if first:
                    h = read_fits_header(first)
                    ra_raw = h.get("OBJCTRA") or h.get("RA") or h.get("CRVAL1", 0)
                    dec_raw = h.get("OBJCTDEC") or h.get("DEC") or h.get("CRVAL2", 0)
                    # Handle sexagesimal strings (e.g. "19 25 29.00") or numeric
                    if isinstance(ra_raw, str) and (' ' in ra_raw or ':' in ra_raw):
                        ra = _hms_to_deg(ra_raw)
                    else:
                        ra = float(ra_raw or 0)
                    if isinstance(dec_raw, str) and (' ' in dec_raw or ':' in dec_raw):
                        dec = _dms_to_deg(dec_raw)
                    else:
                        dec = float(dec_raw or 0)
                    obj      = str(h.get("OBJECT", "?"))
                    n_frames = len(list_fits(lights))

        self.field_ra, self.field_dec = ra, dec
        info = f"{obj}  ·  RA {ra:.4f}°  Dec {dec:+.4f}°  ·  {n_frames} stacked frames"
        self.after(0, lambda: self.lbl_field.configure(text=f"Field: {info}"))

        stars: list[dict] = []
        csv_path = session / "starsv.csv"
        if csv_path.exists():
            with open(csv_path, newline="", encoding="utf-8") as f:
                for row in csv.DictReader(f):
                    try:
                        sr = float(row.get("ra", 0))
                        sd = float(row.get("dec", 0))
                        dist = 0.0
                        if ra:
                            dist = round(
                                (((sr - ra) * 0.7071) ** 2 + (sd - dec) ** 2) ** 0.5, 3
                            )
                        stars.append({
                            "name":     row.get("name", "?"),
                            "ra":       sr, "dec": sd,
                            "mag":      float(row.get("mag", 99) or 99),
                            "var_type": "?", "period": "—", "dist": dist,
                        })
                    except (ValueError, KeyError):
                        continue
            self._log(f"{len(stars)} stars loaded from starsv.csv")
        else:
            self._log("starsv.csv not found — click '↻ Query VSX catalog'")

        self.all_stars = stars
        self.after(0, self._refresh_table)

        proc = session / "process"
        has_registered = any(
            (proc / f"{seq}.seq").exists()
            for seq in ("r_light_", "r_pp_light_")
        )
        hdr: dict = {}
        has_platesolved = False
        if has_registered:
            for seq in ("r_light_", "r_pp_light_"):
                fits_list = sorted(proc.glob(f"{seq.rstrip('_')}*.fit"))
                if fits_list:
                    try:
                        hdr = read_fits_header(fits_list[0])
                        has_platesolved = "CRVAL1" in hdr
                    except Exception:
                        pass
                    break

        has_lc = (proc / "light_curve.dat").exists()
        if has_platesolved and has_lc:
            suggested = "Photometry only (step 5)"
        elif has_platesolved and not has_lc:
            suggested = "Apply reg + plate solve (steps 3–5)"
        elif has_registered:
            suggested = "Plate solve only (step 4–5)"
        else:
            suggested = "Full pipeline (steps 1–5)"

        self.after(0, lambda s=suggested: self.start_from_var.set(s))

        # Best available WCS for in-frame star filtering (plate-solved frame preferred)
        if hdr.get("CRVAL1") and hdr.get("NAXIS1"):
            self.field_wcs_hdr = hdr
            self.field_naxis1 = int(float(hdr.get("NAXIS1", 0)))
            self.field_naxis2 = int(float(hdr.get("NAXIS2", 0)))
        elif h.get("CRVAL1") and h.get("NAXIS1"):
            self.field_wcs_hdr = h
            self.field_naxis1 = int(float(h.get("NAXIS1", 0)))
            self.field_naxis2 = int(float(h.get("NAXIS2", 0)))

    # ── Table ─────────────────────────────────────────────────────────────────

    def _refresh_table(self):
        f = self.filter_var.get().lower()
        try:
            mag_max = float(self.mag_max_var.get())
        except (ValueError, AttributeError):
            mag_max = 99.0
        self.tree.delete(*self.tree.get_children())
        shown = 0
        for s in self.all_stars:
            if f and f not in s["name"].lower() and f not in s["var_type"].lower():
                continue
            if s.get("mag", 99.0) > mag_max:
                continue
            tag = "known" if s.get("var_type", "?") not in ("?", "", " ") else ""
            self.tree.insert("", "end", iid=str(id(s)), tags=(tag,), values=(
                s["name"],
                s.get("var_type", "?"),
                s.get("period", "—"),
                f"{s['mag']:.2f}",
                f"{s.get('dist', 0):.2f}",
            ))
            shown += 1
        self.tree.tag_configure("known", foreground=CYAN)
        self.lbl_count.configure(text=f"{shown} stars")

    def _sort(self, col: str):
        key_map = {"name": "name", "type": "var_type",
                   "period": "period", "mag": "mag", "dist": "dist"}
        key = key_map.get(col, col)
        rev = getattr(self, f"_rev_{col}", False)
        try:
            self.all_stars.sort(
                key=lambda s: float(str(s.get(key, 0)).split()[0] or 0), reverse=rev
            )
        except (ValueError, TypeError):
            self.all_stars.sort(key=lambda s: str(s.get(key, "")), reverse=rev)
        setattr(self, f"_rev_{col}", not rev)
        self._refresh_table()

    def _on_select(self, _=None):
        sel = self.tree.selection()
        if not sel:
            return
        star = next((s for s in self.all_stars if str(id(s)) == sel[0]), None)
        if not star:
            return
        self.selected_star = star
        self.comp_stars = []
        self.lbl_name.configure(text=star["name"])
        self.lbl_type.configure(
            text=f"Type: {star.get('var_type','?')}   Period: {star.get('period','—')}"
        )
        self.lbl_coord.configure(
            text=f"RA {star['ra']:.5f}°  /  Dec {star['dec']:+.5f}°  ·  mag {star['mag']:.2f}"
        )
        self.run_btn.configure(state="normal", bg="#3a7bd5", fg=FG)

    # ── VSX catalog ───────────────────────────────────────────────────────────

    def _fetch_vsx(self):
        if not self.field_ra:
            messagebox.showinfo("Info", "Load a session first.")
            return
        if not HAS_REQUESTS:
            messagebox.showerror("Missing",
                                 "requests not available.\n"
                                 "Run: test_connections.py for diagnostics.")
            return
        self._log("Querying VizieR VSX catalog…")
        self._btn_state(False)
        threading.Thread(target=self._fetch_vsx_bg, daemon=True).start()

    def _fetch_vsx_bg(self):
        radius = max(SEESTAR["fov_w"], SEESTAR["fov_h"]) * 0.75
        stars = query_vsx(self.field_ra, self.field_dec, radius)
        if not stars:
            self.after(0, lambda: (
                self._log("VSX: no results — check connection with test_connections.py"),
                self._btn_state(True),
            ))
            return
        for s in stars:
            s["dist"] = round(
                (((s["ra"] - self.field_ra) * 0.7071) ** 2
                 + (s["dec"] - self.field_dec) ** 2) ** 0.5, 3
            )

        # Keep only stars inside the alt-az safe zone: the inscribed circle of
        # the frame (radius = min(W,H)/2).  With alt-az tracking, field rotation
        # means only this circle is guaranteed to be covered by ALL registered
        # frames — stars outside it can land in the black rotation corners.
        if self.field_wcs_hdr and self.field_naxis1 and self.field_naxis2:
            filtered = stars_in_safe_circle(stars, self.field_wcs_hdr,
                                            self.field_naxis1, self.field_naxis2,
                                            margin=50)
            safe_r = min(self.field_naxis1, self.field_naxis2) / 2 - 50
            self._log(f"VSX: {len(stars)} in search area → "
                      f"{len(filtered)} inside alt-az safe zone "
                      f"(inscribed circle r={safe_r:.0f} px)")
            stars = filtered
        else:
            self._log(f"VSX: {len(stars)} found (no WCS yet — plate-solve for exact filtering)")

        # Persist filtered list so it survives restart
        if self.session_dir and stars:
            csv_path = self.session_dir / "starsv.csv"
            with open(csv_path, "w", encoding="utf-8", newline="") as f:
                w = csv.writer(f)
                w.writerow(["name", "ra", "dec", "mag", "var_type", "period"])
                for s in stars:
                    w.writerow([s["name"], s["ra"], s["dec"],
                                 s.get("mag", 99), s.get("var_type", "?"),
                                 s.get("period", "—")])

        # Write in-frame stars to Siril's user DSO catalogue (for manual Annotate in GUI)
        cat_path = find_siril_user_catalogue()
        if cat_path and stars:
            n_new = update_siril_catalogue(stars, cat_path)
            if n_new:
                self._log(f"Siril catalogue: {n_new} new star(s) added to {cat_path.name}")
                self._log("  → in Siril: Outils > Astrométrie > Annoter to see them")

        self.all_stars = sorted(stars, key=lambda s: s["dist"])
        self.after(0, lambda: (
            self._refresh_table(),
            self._log(f"VSX: {len(stars)} variable stars in field"),
            self._btn_state(True),
        ))

    # ── Pipeline ──────────────────────────────────────────────────────────────

    def _run_pipeline(self):
        if not self.selected_star:
            return
        self.run_btn.configure(state="disabled", bg=SURFACE, fg=FG2)
        self.prog_var.set(0)
        self.log.delete("1.0", "end")
        threading.Thread(target=self._pipeline_bg, daemon=True).start()

    def _get_start_step(self) -> int:
        label = self.start_from_var.get()
        if "step 5" in label or "Photometry" in label:
            return 5
        if "step 4" in label or "Plate solve" in label:
            return 4
        if "step 3" in label or "Apply reg" in label:
            return 3
        return 1

    def _pipeline_bg(self):
        filt_key = self.obs_filter_var.get()
        filt_code, filt_note = FILTER_OPTIONS.get(filt_key, ("CV", ""))
        config = {
            "session":    self.session_dir,
            "star":       self.selected_star,
            "comp_stars": list(self.comp_stars),
            "nstars":     self.nstars_var.get(),
            "start_step": self._get_start_step(),
            "dark_dir":   Path(self.dark_var.get()) if self.dark_var.get().strip() else None,
            "flat_dir":   Path(self.flat_var.get()) if self.flat_var.get().strip() else None,
            "bias_dir":   Path(self.bias_var.get()) if self.bias_var.get().strip() else None,
            "filt_code":  filt_code,
            "filt_note":  filt_note,
            "runner":     self.runner,
            "phot_aperture":  self.ap_var.get(),
            "phot_inner":     self.inner_var.get(),
            "phot_outer":     self.outer_var.get(),
            "phot_dyn_ratio": self.dyn_var.get(),
            "phot_max_val":   self.maxval_var.get(),
        }
        run_pipeline(config, self._log, self._prog, self._done)

    # ── Thread-safe helpers ───────────────────────────────────────────────────

    def _log(self, msg: str):
        def _do():
            self.log.insert("end", msg + "\n")
            self.log.see("end")
        self.after(0, _do)

    def _prog(self, pct: int, label: str):
        self.after(0, lambda: (
            self.prog_var.set(pct),
            self.lbl_prog.configure(text=label, fg=CYAN),
        ))

    def _btn_state(self, enabled: bool):
        if self.selected_star:
            self.after(0, lambda: self.run_btn.configure(
                state="normal" if enabled else "disabled",
                bg=BLUE if enabled else SURFACE,
                fg=FG if enabled else FG2,
            ))

    def _done(self, success: bool, msg: str):
        self.after(0, lambda: self._done_ui(success, msg))

    def _done_ui(self, ok: bool, msg: str):
        self.run_btn.configure(state="normal", bg=BLUE, fg=FG)
        self.prog_var.set(100 if ok else 0)
        if ok:
            self.lbl_prog.configure(text="Done!", fg=GREEN)
            self._log(f"✓  {msg}")
            messagebox.showinfo("Done!",
                                f"Results saved in:\n{msg}\n\n"
                                "Files: light_curve.dat  ·  aavso.csv\n"
                                "       photometry.csv  ·  light_curve.png\n"
                                "       pipeline.log")
        else:
            self.lbl_prog.configure(text=f"Failed: {msg}", fg=RED)
            self._log(f"✗  {msg}")
            messagebox.showerror("Pipeline error", msg)

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--session", help="Session folder (auto-detected if omitted)")
    args, _ = p.parse_known_args()
    session = Path(args.session) if args.session else None
    App(session_dir=session).mainloop()
