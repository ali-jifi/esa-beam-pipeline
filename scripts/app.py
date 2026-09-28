# tkinter front end for the beam pipeline, spectrogram, and candidate cli scripts
import ast
import calendar
import csv
import json
import math
import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
import webbrowser
from datetime import date, datetime, timedelta
from pathlib import Path
from tkinter import filedialog, messagebox, simpledialog, ttk
from tkinter.scrolledtext import ScrolledText

from PIL import Image, ImageTk

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
FIGURES = ROOT / "figures"
CANDIDATES = ROOT / "candidates"
PRESETS_FILE = ROOT / "configs" / "ui_presets.json"
PROBES = "abcde"
FIG_NAMES = ("overview", "snapshots", "histograms", "threshold_compare")
INT_KEYS = {"min_consecutive", "coherent_min_bins"}

# key, cli flag, label, tip, basic
PARAMS = [
    ("score_threshold", "--score-threshold", "score threshold", "beam score threshold", True),
    ("asymmetry_min", "--asym-threshold", "asymmetry", "asymmetry score-ramp anchor", False),
    ("width_max", "--width-threshold", "width", "width score-ramp anchor", False),
    ("para_to_omni_min", "--p2o-threshold", "para/omni", "para-to-omni score-ramp anchor", False),
    ("min_coverage", "--min-coverage", "min coverage", "min pa cone solid-angle coverage", False),
    ("n_sigma_lo", "--n-sigma-lo", "n sigma lo", "lo poisson bar, runs form here, needs a neighbor", False),
    ("n_sigma_hi", "--n-sigma-hi", "n sigma hi", "hi poisson bar, isolated runs stand alone", False),
    ("beam_e_max", "--beam-e-max", "beam E max (eV)", "band-energy ceiling for any accepted run", False),
    ("coherent_asym_min", "--coherent-asym-min", "coherent |asym| min", "per-bin |asym| threshold for coherent run", False),
    ("coherent_dir_min", "--coherent-dir-min", "coherent R min", "per-bin dominant-cone/perp (R) threshold", False),
    ("coherent_min_bins", "--coherent-min-bins", "coherent min bins", "min adjacent bins for coherent beam", False),
    ("peak_prom_min", "--peak-prom-min", "peak prom min", "log10 prominence for spectral line, 0.3 = 2x", False),
    ("peak_width_max", "--peak-width-max", "peak width max", "max line fwhm in bins", False),
]
# key, cli flag, label, tip, default
EXTRA = [
    ("energy_cutoff", "--energy-cutoff", "energy cutoff (eV)", "low-energy cutoff", 30.0),
    ("min_consecutive", "--min-consecutive", "min consecutive", "min consecutive beam steps to keep, 1 = keep isolated", 1),
]


def classifier_defaults() -> dict:
    # field defaults straight from the dataclass source, skips the 8s pyspedas import
    tree = ast.parse((SRC / "esa_plotting" / "beam_pipeline.py").read_text(encoding="utf-8-sig"))
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "ClassifierParams")
    return {n.target.id: ast.literal_eval(n.value) for n in cls.body
            if isinstance(n, ast.AnnAssign) and n.value is not None}


def builtin_presets(d: dict) -> dict:
    strict = {k: d[k] for k, *_ in PARAMS} | {k: v for k, _, _, _, v in EXTRA}
    # mirrors make_profile in generate_candidates
    survey = strict | {"n_sigma_lo": 0.5 * d["n_sigma_lo"], "n_sigma_hi": 0.5 * d["n_sigma_hi"],
                       "coherent_dir_min": 0.5 * d["coherent_dir_min"], "coherent_min_bins": 2}
    return {"strict (default)": strict, "survey": survey}


def load_user_presets() -> dict:
    try:
        return json.loads(PRESETS_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def load_events() -> dict:
    out = {}
    for f in sorted(ROOT.glob("events_*.csv")):
        with open(f, newline="") as fh:
            out[f.stem.removeprefix("events_")] = list(csv.DictReader(fh))
    return out


def clean(key: str, s: str) -> str:
    return str(int(float(s))) if key in INT_KEYS else str(float(s))


def flat(rec: dict) -> dict:
    out = {}
    for k, v in rec.items():
        if isinstance(v, dict):
            out.update({f"{k}.{kk}": vv for kk, vv in v.items()})
        else:
            out[k] = v
    return out


def open_path(p: Path):
    if hasattr(os, "startfile"):
        os.startfile(p)
    else:
        webbrowser.open(p.as_uri())


class Tip:
    # hover tip, delayed show, hides on leave or click, polls the pointer while up
    def __init__(self, w, text):
        self.w, self.text, self.top, self.job = w, text, None, None
        w.bind("<Enter>", self._arm, add="+")
        for ev in ("<Leave>", "<ButtonPress>", "<Destroy>"):
            w.bind(ev, self.hide, add="+")

    def _arm(self, _):
        self.hide()
        self.job = self.w.after(500, self.show)

    def show(self):
        self.job = None
        if not self._inside():
            return
        self.top = tk.Toplevel(self.w)
        self.top.wm_overrideredirect(True)
        self.top.wm_geometry(f"+{self.w.winfo_rootx()}+{self.w.winfo_rooty() + self.w.winfo_height() + 4}")
        tk.Label(self.top, text=self.text, bg="#ffffe0", relief="solid", bd=1, padx=4, pady=2).pack()
        self.job = self.w.after(300, self._watch)

    def _watch(self):
        # combobox popups and dialogs grab the pointer and swallow the leave event
        self.job = None
        if self._inside():
            self.job = self.w.after(300, self._watch)
        else:
            self.hide()

    def _inside(self):
        try:
            x, y = self.w.winfo_pointerxy()
            x0, y0 = self.w.winfo_rootx(), self.w.winfo_rooty()
            return x0 <= x < x0 + self.w.winfo_width() and y0 <= y < y0 + self.w.winfo_height()
        except tk.TclError:
            return False

    def hide(self, *_):
        if self.job:
            self.w.after_cancel(self.job)
            self.job = None
        if self.top:
            self.top.destroy()
            self.top = None


class CalendarPopup(tk.Toplevel):
    def __init__(self, master, var: tk.StringVar):
        super().__init__(master)
        self.var = var
        self.title("pick date")
        self.resizable(False, False)
        self.transient(master.winfo_toplevel())
        try:
            d = datetime.strptime(var.get(), "%Y-%m-%d").date()
        except ValueError:
            d = date.today()
        self.y, self.m = d.year, d.month
        hdr = ttk.Frame(self)
        hdr.pack(fill="x", padx=4, pady=4)
        ttk.Button(hdr, text="<", width=3, command=lambda: self.shift(-1)).pack(side="left")
        ttk.Button(hdr, text=">", width=3, command=lambda: self.shift(1)).pack(side="right")
        self.lbl = ttk.Label(hdr, anchor="center")
        self.lbl.pack(fill="x")
        self.days = ttk.Frame(self)
        self.days.pack(padx=4, pady=4)
        self.draw()

    def shift(self, n):
        self.y, m0 = divmod(self.y * 12 + self.m - 1 + n, 12)
        self.m = m0 + 1
        self.draw()

    def draw(self):
        for w in self.days.winfo_children():
            w.destroy()
        self.lbl.config(text=f"{calendar.month_name[self.m]} {self.y}")
        for c, name in enumerate(calendar.day_abbr):
            ttk.Label(self.days, text=name[:2], width=4, anchor="center").grid(row=0, column=c)
        for r, week in enumerate(calendar.monthcalendar(self.y, self.m), start=1):
            for c, day in enumerate(week):
                if day:
                    ttk.Button(self.days, text=str(day), width=4,
                               command=lambda d=day: self.pick(d)).grid(row=r, column=c)

    def pick(self, d):
        self.var.set(f"{self.y:04d}-{self.m:02d}-{d:02d}")
        self.destroy()


def row(parent, r, label, widget, tip=None, col=0):
    ttk.Label(parent, text=label).grid(row=r, column=2 * col, sticky="w", padx=4, pady=2)
    widget.grid(row=r, column=2 * col + 1, sticky="w", padx=4, pady=2)
    if tip:
        Tip(widget, tip)


def date_field(parent, var):
    f = ttk.Frame(parent)
    ttk.Entry(f, textvariable=var, width=11).pack(side="left")
    ttk.Button(f, text="...", width=3, command=lambda: CalendarPopup(parent, var)).pack(side="left")
    return f


class IntervalPicker(ttk.LabelFrame):
    # probe + start/end date and time, or a known interval from an events csv
    def __init__(self, master, events: dict):
        super().__init__(master, text="interval")
        self.events = events
        self.rows = []
        self.probe = tk.StringVar(value="a")
        self.d0, self.t0 = tk.StringVar(value="2019-05-01"), tk.StringVar(value="00:00")
        self.d1, self.t1 = tk.StringVar(value="2019-05-02"), tk.StringVar(value="00:00")
        self.d0.trace_add("write", self._roll_end)
        row(self, 0, "probe", ttk.Combobox(self, textvariable=self.probe, values=list(PROBES),
                                          width=4, state="readonly"), "themis probe a-e")
        for r, (label, d, t) in enumerate((("start", self.d0, self.t0), ("end", self.d1, self.t1)), start=1):
            f = ttk.Frame(self)
            date_field(f, d).pack(side="left")
            ttk.Entry(f, textvariable=t, width=6).pack(side="left", padx=(6, 0))
            row(self, r, label, f, "date and ut time HH:MM, 00:00 to 00:00 next day = full day")
        self.src = ttk.Combobox(self, values=list(events), width=16, state="readonly")
        self.src.bind("<<ComboboxSelected>>", self._fill_list)
        row(self, 3, "known intervals", self.src, "pick a csv from the repo events lists")
        self.evt = ttk.Combobox(self, width=44, state="readonly")
        self.evt.bind("<<ComboboxSelected>>", self._pick)
        self.evt.grid(row=4, column=0, columnspan=2, sticky="w", padx=4, pady=2)

    def _roll_end(self, *_):
        try:
            a = datetime.strptime(self.d0.get(), "%Y-%m-%d")
            b = datetime.strptime(self.d1.get(), "%Y-%m-%d")
        except ValueError:
            return
        if b < a:
            self.d1.set((a + timedelta(days=1)).strftime("%Y-%m-%d"))

    def _fill_list(self, _):
        self.rows = self.events[self.src.get()]
        self.evt["values"] = [f"{r['probe']}  {r['trange_start']} -> {r['trange_end']}  {r['interval_tag']}"
                              for r in self.rows]
        self.evt.set("")

    def _pick(self, _):
        r = self.rows[self.evt.current()]
        self.probe.set(r["probe"])
        for s, d, t in ((r["trange_start"], self.d0, self.t0), (r["trange_end"], self.d1, self.t1)):
            day, _, tm = s.partition("/")
            d.set(day)
            t.set(tm or "00:00")

    def trange(self) -> list[str]:
        out, stamps = [], []
        for d, t in ((self.d0, self.t0), (self.d1, self.t1)):
            day, tm = d.get().strip(), t.get().strip()
            stamps.append(datetime.strptime(f"{day} {tm}", "%Y-%m-%d %H:%M"))
            out.append(day if tm == "00:00" else f"{day}/{tm}")
        if stamps[1] <= stamps[0]:
            raise ValueError("end must be after start")
        return out


class Runner:
    # runs a script in a subprocess, streams its output into the log
    def __init__(self, root, log: ScrolledText, status: tk.StringVar):
        self.root, self.log, self.status = root, log, status
        self.proc, self.q, self.on_done, self.cancelled = None, queue.Queue(), None, False

    @property
    def busy(self):
        return self.proc is not None and self.proc.poll() is None

    def start(self, cmd, env, on_done):
        if self.busy:
            messagebox.showwarning("busy", "a run is already going")
            return
        cmd = [str(c) for c in cmd]
        self.log.delete("1.0", "end")
        self.append("$ " + " ".join(cmd) + "\n\n")
        self.on_done, self.cancelled = on_done, False
        self.proc = subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                                     errors="replace", bufsize=1,
                                     creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        threading.Thread(target=self._pump, daemon=True).start()
        self.status.set("running")
        self.root.after(100, self._poll)

    def _pump(self):
        for line in self.proc.stdout:
            self.q.put(line)
        self.q.put(None)

    def _poll(self):
        while True:
            try:
                line = self.q.get_nowait()
            except queue.Empty:
                break
            if line is None:
                rc = self.proc.wait()
                self.status.set("cancelled" if self.cancelled else "done" if rc == 0 else f"failed, exit {rc}")
                self.append(f"\n[exit {rc}]\n")
                self.on_done(rc)
                return
            self.append(line)
        self.root.after(100, self._poll)

    def append(self, s):
        self.log.insert("end", s)
        self.log.see("end")

    def cancel(self):
        if self.busy:
            self.cancelled = True
            self.proc.terminate()


class ImageView(tk.Canvas):
    # fits by default, wheel zooms at the cursor, drag pans, double click refits
    def __init__(self, master):
        super().__init__(master, highlightthickness=0,
                         bg=ttk.Style().lookup("TFrame", "background") or "#f0f0f0")
        self.img, self.photo, self.scale, self.fit_scale, self.ox, self.oy = None, None, 1.0, 1.0, 0.0, 0.0
        self.fitted, self._job, self._from = True, None, None
        self.bind("<Configure>", self._resized)
        self.bind("<MouseWheel>", lambda e: self.zoom(1.25 if e.delta > 0 else 0.8, e.x, e.y))
        self.bind("<Button-4>", lambda e: self.zoom(1.25, e.x, e.y))
        self.bind("<Button-5>", lambda e: self.zoom(0.8, e.x, e.y))
        self.bind("<ButtonPress-1>", self._press)
        self.bind("<B1-Motion>", self._drag)
        self.bind("<Double-Button-1>", lambda e: self.fit())

    def set_image(self, img):
        self.img = img
        self.fit()

    def _size(self):
        return max(self.winfo_width(), 50), max(self.winfo_height(), 50)

    def fit(self):
        self.fitted = True
        if self.img is None:
            self.delete("all")
            return
        w, h = self._size()
        self.scale = min(w / self.img.width, h / self.img.height, 1.0)
        self.redraw()

    def zoom(self, f, x, y):
        if self.img is None:
            return
        new = min(max(self.scale * f, self.fit_scale), 8.0)
        # keep the image point under the cursor where it is
        self.ox = x - (x - self.ox) * new / self.scale
        self.oy = y - (y - self.oy) * new / self.scale
        self.scale, self.fitted = new, new <= self.fit_scale
        self.redraw()

    def _press(self, e):
        self._from = (e.x, e.y)

    def _drag(self, e):
        if self._from and self.img is not None:
            self.ox += e.x - self._from[0]
            self.oy += e.y - self._from[1]
            self._from = (e.x, e.y)
            self.redraw()

    def _resized(self, _):
        if self._job:
            self.after_cancel(self._job)
        self._job = self.after(80, self.fit if self.fitted else self.redraw)

    def redraw(self):
        self._job = None
        if self.img is None:
            return
        w, h = self._size()
        self.fit_scale = min(w / self.img.width, h / self.img.height, 1.0)
        iw, ih = self.img.width * self.scale, self.img.height * self.scale
        # center when smaller than the canvas, else keep the edges inside it
        self.ox = (w - iw) / 2 if iw <= w else min(0.0, max(w - iw, self.ox))
        self.oy = (h - ih) / 2 if ih <= h else min(0.0, max(h - ih, self.oy))
        # render only the visible part so deep zooms stay cheap
        x0, y0 = int(max(0.0, -self.ox / self.scale)), int(max(0.0, -self.oy / self.scale))
        x1 = min(self.img.width, math.ceil((w - self.ox) / self.scale))
        y1 = min(self.img.height, math.ceil((h - self.oy) / self.scale))
        size = (max(1, round((x1 - x0) * self.scale)), max(1, round((y1 - y0) * self.scale)))
        res = Image.Resampling.LANCZOS if self.scale < 1 else Image.Resampling.BILINEAR
        self.photo = ImageTk.PhotoImage(self.img.crop((x0, y0, x1, y1)).resize(size, res))
        self.delete("all")
        self.create_image(self.ox + x0 * self.scale, self.oy + y0 * self.scale, image=self.photo, anchor="nw")


class ResultsPanel(ttk.Frame):
    # figure viewer plus table/text tabs, loads a run from its file prefix
    def __init__(self, master):
        super().__init__(master)
        self.figs, self.run_paths = {}, []
        top = ttk.Frame(self)
        top.pack(fill="x", pady=2)
        ttk.Label(top, text="previous runs").pack(side="left")
        self.runs = ttk.Combobox(top, state="readonly", width=44)
        self.runs.pack(side="left", padx=4)
        self.runs.bind("<<ComboboxSelected>>", lambda e: self.load_prefix(self.run_paths[self.runs.current()]))
        ttk.Button(top, text="refresh", command=self.refresh_runs).pack(side="left")
        mid = ttk.Frame(self)
        mid.pack(fill="x", pady=2)
        ttk.Label(mid, text="figure").pack(side="left")
        self.fig = ttk.Combobox(mid, state="readonly", width=32)
        self.fig.pack(side="left", padx=4)
        self.fig.bind("<<ComboboxSelected>>", self.show_fig)
        ttk.Button(mid, text="open", command=lambda: self._open(False)).pack(side="left")
        ttk.Button(mid, text="open folder", command=lambda: self._open(True)).pack(side="left", padx=4)
        ttk.Button(mid, text="fit", command=lambda: self.view.fit()).pack(side="left")
        ttk.Label(mid, text="wheel zooms, drag pans, double click refits").pack(side="right", padx=4)
        pane = ttk.PanedWindow(self, orient="vertical")
        pane.pack(fill="both", expand=True)
        # ttk panes size children by request so the empty canvas starts tiny, give it 70% once
        pane.bind("<Configure>", lambda e: (pane.sashpos(0, int(e.height * 0.7)), pane.unbind("<Configure>"))
                  if e.height > 100 else None)
        self.view = ImageView(pane)
        pane.add(self.view, weight=3)
        tabs = ttk.Notebook(pane)
        pane.add(tabs, weight=1)
        tf = ttk.Frame(tabs)
        tabs.add(tf, text="table")
        self.tree = ttk.Treeview(tf, show="headings")
        ys = ttk.Scrollbar(tf, orient="vertical", command=self.tree.yview)
        xs = ttk.Scrollbar(tf, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        ys.pack(side="right", fill="y")
        xs.pack(side="bottom", fill="x")
        self.tree.pack(fill="both", expand=True)
        self.text = ScrolledText(tabs, font="TkFixedFont", height=8)
        tabs.add(self.text, text="diagnose")

    def refresh_runs(self):
        pngs = sorted(FIGURES.glob("*/*_overview.png")) + sorted(FIGURES.glob("*.png"))
        self.run_paths = [p.with_name(p.name.removesuffix("_overview.png")) if p.name.endswith("_overview.png")
                          else p for p in pngs]
        self.runs["values"] = [str(p.relative_to(FIGURES)) for p in self.run_paths]

    def load_prefix(self, prefix: Path):
        # pipeline prefix like figures/2019-05-01/tha_beam_2019-05-01, or a bare png
        if prefix.suffix == ".png":
            self.show_figs({prefix.stem: prefix} if prefix.exists() else {})
            self.set_table([])
            self.set_text("")
            return
        figs = {n: p for n in FIG_NAMES if (p := prefix.with_name(f"{prefix.name}_{n}.png")).exists()}
        self.show_figs(figs)
        beams = prefix.with_name(f"{prefix.name}_beams.csv")
        if beams.exists():
            with open(beams, newline="") as fh:
                self.set_table(list(csv.DictReader(fh)))
        else:
            self.set_table([])
        diag = prefix.with_name(f"{prefix.name}_diagnose.txt")
        self.set_text(diag.read_text(encoding="utf-8", errors="replace") if diag.exists() else "")

    def load_candidates(self, jsonl: Path):
        if not jsonl.exists():
            self.show_figs({})
            self.set_table([])
            return
        recs = [json.loads(l) for l in jsonl.read_text(encoding="utf-8").splitlines() if l.strip()]
        self.set_table([flat(r) for r in recs])
        cut = jsonl.parent / "cutouts"
        self.show_figs({r["candidate_id"]: p for r in recs if (p := cut / f"{r['candidate_id']}.png").exists()})
        self.set_text("")

    def show_figs(self, figs: dict):
        self.figs = figs
        self.fig["values"] = list(figs)
        if figs:
            self.fig.current(0)
            self.show_fig()
        else:
            self.fig.set("")
            self.view.set_image(None)

    def show_fig(self, *_):
        p = self.figs.get(self.fig.get())
        if p:
            with Image.open(p) as im:
                self.view.set_image(im.copy())

    def _open(self, folder: bool):
        p = self.figs.get(self.fig.get())
        if p:
            open_path(p.parent if folder else p)

    def set_table(self, rows: list):
        self.tree.delete(*self.tree.get_children())
        cols = list(dict.fromkeys(k for r in rows for k in r))
        self.tree["columns"] = cols
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=max(70, 8 * len(c)), stretch=False)
        for r in rows:
            self.tree.insert("", "end", values=[f"{v:.4g}" if isinstance(v, float) else v
                                                for v in (r.get(c, "") for c in cols)])

    def set_text(self, s: str):
        self.text.delete("1.0", "end")
        self.text.insert("1.0", s)


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("esa_plotting")
        self.geometry("1500x950")
        self.defaults = classifier_defaults()
        self.events = load_events()
        self.presets = builtin_presets(self.defaults) | load_user_presets()
        self.data_dir = tk.StringVar(value=os.environ.get("THM_DATA_DIR") or str(ROOT / "data"))
        self.status = tk.StringVar(value="idle")

        top = ttk.Frame(self, padding=4)
        top.pack(fill="x")
        ttk.Label(top, text="data dir").pack(side="left")
        e = ttk.Entry(top, textvariable=self.data_dir, width=70)
        e.pack(side="left", padx=4)
        Tip(e, "where themis cdfs live, missing files download from the themis server on first run")
        ttk.Button(top, text="browse", command=self.pick_data_dir).pack(side="left")
        ttk.Label(top, textvariable=self.status).pack(side="right", padx=8)

        vpane = ttk.PanedWindow(self, orient="vertical")
        vpane.pack(fill="both", expand=True)
        hpane = ttk.PanedWindow(vpane, orient="horizontal")
        vpane.add(hpane, weight=4)
        self.tabs = ttk.Notebook(hpane)
        hpane.add(self.tabs, weight=0)
        self.results = ResultsPanel(hpane)
        hpane.add(self.results, weight=1)
        self.log = ScrolledText(vpane, height=10, font="TkFixedFont")
        vpane.add(self.log, weight=1)
        self.runner = Runner(self, self.log, self.status)

        self.build_pipeline_tab()
        self.build_spectro_tab()
        self.build_cand_tab()
        self.results.refresh_runs()
        self.protocol("WM_DELETE_WINDOW", self.close)

    def env(self) -> dict:
        e = os.environ.copy()
        e["THM_DATA_DIR"] = self.data_dir.get().strip()
        e["PYTHONPATH"] = str(SRC) + os.pathsep + e.get("PYTHONPATH", "")
        e["PYTHONIOENCODING"] = "utf-8"
        e["MPLBACKEND"] = "Agg"
        return e

    def pick_data_dir(self):
        d = filedialog.askdirectory(initialdir=self.data_dir.get())
        if d:
            self.data_dir.set(d)

    def run_buttons(self, parent, cmd):
        bf = ttk.Frame(parent)
        bf.pack(fill="x", pady=6)
        ttk.Button(bf, text="run", command=cmd).pack(side="left")
        ttk.Button(bf, text="cancel", command=self.runner.cancel).pack(side="left", padx=4)

    def build_pipeline_tab(self):
        f = ttk.Frame(self.tabs, padding=6)
        self.tabs.add(f, text="beam pipeline")
        self.pp_int = IntervalPicker(f, self.events)
        self.pp_int.pack(fill="x")

        pf = ttk.Frame(f)
        pf.pack(fill="x", pady=4)
        ttk.Label(pf, text="preset").pack(side="left")
        self.preset = ttk.Combobox(pf, values=list(self.presets), state="readonly", width=18)
        self.preset.pack(side="left", padx=4)
        self.preset.bind("<<ComboboxSelected>>", lambda e: self.apply_preset(self.preset.get()))
        ttk.Button(pf, text="save as", command=self.save_preset).pack(side="left")

        self.pvars = {k: tk.StringVar() for k, *_ in EXTRA + PARAMS}
        basic = ttk.LabelFrame(f, text="parameters")
        basic.pack(fill="x")
        for r, (key, _, label, tip, _) in enumerate(EXTRA + [p for p in PARAMS if p[4]]):
            row(basic, r, label, ttk.Entry(basic, textvariable=self.pvars[key], width=8), tip)

        self.show_adv = tk.BooleanVar()
        ttk.Checkbutton(f, text="show advanced", variable=self.show_adv, command=self.toggle_adv).pack(anchor="w", pady=2)
        self.adv = ttk.LabelFrame(f, text="advanced")
        for i, (key, _, label, tip, _) in enumerate([p for p in PARAMS if not p[4]]):
            row(self.adv, i // 2, label, ttk.Entry(self.adv, textvariable=self.pvars[key], width=8), tip, col=i % 2)

        self.opt = ttk.LabelFrame(f, text="options")
        self.opt.pack(fill="x", pady=4)
        self.tcv = tk.StringVar(value="1.0 1.2 1.5 2.0")
        row(self.opt, 0, "threshold compare R", ttk.Entry(self.opt, textvariable=self.tcv, width=20),
            "coherent R values for the threshold-comparison plot, space separated")
        self.no_plots = tk.BooleanVar()
        ttk.Checkbutton(self.opt, text="skip plots", variable=self.no_plots).grid(row=1, column=0, columnspan=2, sticky="w", padx=4)
        self.diag = tk.BooleanVar()
        ttk.Checkbutton(self.opt, text="diagnose dump", variable=self.diag).grid(row=2, column=0, columnspan=2, sticky="w", padx=4)
        self.diag0, self.diag1 = tk.StringVar(), tk.StringVar()
        df = ttk.Frame(self.opt)
        ttk.Entry(df, textvariable=self.diag0, width=6).pack(side="left")
        ttk.Label(df, text="to").pack(side="left", padx=4)
        ttk.Entry(df, textvariable=self.diag1, width=6).pack(side="left")
        row(self.opt, 3, "diagnose UT window", df, "HH:MM HH:MM, blank = whole run window")

        self.run_buttons(f, self.run_pipeline)
        self.preset.current(0)
        self.apply_preset(self.preset.get())

    def toggle_adv(self):
        if self.show_adv.get():
            self.adv.pack(fill="x", before=self.opt)
        else:
            self.adv.pack_forget()

    def apply_preset(self, name):
        for k, v in self.presets[name].items():
            if k in self.pvars:
                self.pvars[k].set(str(v))

    def save_preset(self):
        name = simpledialog.askstring("save preset", "preset name", parent=self)
        if not name:
            return
        try:
            vals = {k: int(float(v.get())) if k in INT_KEYS else float(v.get()) for k, v in self.pvars.items()}
        except ValueError as e:
            messagebox.showerror("bad input", str(e))
            return
        user = load_user_presets() | {name: vals}
        PRESETS_FILE.parent.mkdir(exist_ok=True)
        PRESETS_FILE.write_text(json.dumps(user, indent=2), encoding="utf-8")
        self.presets[name] = vals
        self.preset["values"] = list(self.presets)
        self.preset.set(name)

    def run_pipeline(self):
        try:
            trange = self.pp_int.trange()
            args = [(flag, clean(key, self.pvars[key].get())) for key, flag, *_ in EXTRA + PARAMS]
            tcv = [str(float(x)) for x in self.tcv.get().split()]
        except ValueError as e:
            messagebox.showerror("bad input", str(e))
            return
        probe = self.pp_int.probe.get()
        cmd = [sys.executable, "-u", SCRIPTS / "run_beam_pipeline.py", "--probe", probe, "--trange", *trange]
        for flag, v in args:
            cmd += [flag, v]
        if tcv:
            cmd += ["--threshold-compare-values", *tcv]
        if self.no_plots.get():
            cmd.append("--no-plots")
        if self.diag.get():
            cmd.append("--diagnose")
            win = [self.diag0.get().strip(), self.diag1.get().strip()]
            if all(win):
                cmd += win
        day = trange[0].split("/")[0]
        tag = "-".join(s.split("/")[1].replace(":", "") for s in trange if "/" in s)
        prefix = FIGURES / day / f"th{probe}_beam_{day}{'_' + tag if tag else ''}"
        self.runner.start(cmd, self.env(), lambda rc: (self.results.load_prefix(prefix), self.results.refresh_runs()))

    def build_spectro_tab(self):
        f = ttk.Frame(self.tabs, padding=6)
        self.tabs.add(f, text="spectrogram")
        g = ttk.LabelFrame(f, text="ion eflux")
        g.pack(fill="x")
        self.sp_date = tk.StringVar(value="2019-05-01")
        row(g, 0, "date", date_field(g, self.sp_date), "one full day")
        self.sp_probes = {p: tk.BooleanVar(value=p == "a") for p in PROBES}
        pf = ttk.Frame(g)
        for p in PROBES:
            ttk.Checkbutton(pf, text=p, variable=self.sp_probes[p]).pack(side="left")
        row(g, 1, "probes", pf, "one probe = single panel, several = stacked")
        self.run_buttons(f, self.run_spectro)

    def run_spectro(self):
        d = self.sp_date.get().strip()
        probes = [p for p in PROBES if self.sp_probes[p].get()]
        try:
            datetime.strptime(d, "%Y-%m-%d")
            if not probes:
                raise ValueError("pick at least one probe")
        except ValueError as e:
            messagebox.showerror("bad input", str(e))
            return
        ymd = d.replace("-", "")
        if len(probes) == 1:
            cmd = [sys.executable, "-u", SCRIPTS / "plot_single_probe.py", "--probe", probes[0], "--date", d]
            out = FIGURES / f"th{probes[0]}_peif_eflux_{ymd}.png"
        else:
            cmd = [sys.executable, "-u", SCRIPTS / "plot_multi_probe.py", "--date", d, "--probes", "".join(probes)]
            out = FIGURES / f"multi_probe_eflux_{ymd}.png"
        self.runner.start(cmd, self.env(), lambda rc: (self.results.load_prefix(out), self.results.refresh_runs()))

    def build_cand_tab(self):
        f = ttk.Frame(self.tabs, padding=6)
        self.tabs.add(f, text="candidates")
        self.cd_int = IntervalPicker(f, self.events)
        self.cd_int.pack(fill="x")
        g = ttk.LabelFrame(f, text="generation")
        g.pack(fill="x", pady=4)
        self.cd_tag = tk.StringVar(value="unknown")
        row(g, 0, "interval tag", ttk.Combobox(g, textvariable=self.cd_tag, state="readonly", width=12,
                                               values=["beam_period", "plasma_sheet", "unknown"]),
            "label context stored on every record")
        self.cd_profile = tk.StringVar(value="strict")
        row(g, 1, "profile", ttk.Combobox(g, textvariable=self.cd_profile, state="readonly", width=12,
                                          values=["strict", "survey"]),
            "strict = production params, survey = sigma bars and R gate at 0.5x")
        self.cd_ecut = tk.StringVar(value="30.0")
        row(g, 2, "energy cutoff (eV)", ttk.Entry(g, textvariable=self.cd_ecut, width=8), "low-energy cutoff")
        self.cd_dtype = tk.StringVar(value="peif")
        row(g, 3, "datatype", ttk.Combobox(g, textvariable=self.cd_dtype, state="readonly", width=12,
                                           values=["peif", "peir"]),
            "peif = full survey ~100 s cadence, peir = reduced 3.2 s")
        self.cd_out = tk.StringVar(value=str(CANDIDATES))
        of = ttk.Frame(g)
        ttk.Entry(of, textvariable=self.cd_out, width=30).pack(side="left")
        ttk.Button(of, text="...", width=3, command=self.pick_cand_dir).pack(side="left")
        row(g, 4, "output dir", of, "jsonl + cutouts land here, an existing jsonl is kept and the run is skipped")
        self.cd_nocut = tk.BooleanVar()
        ttk.Checkbutton(g, text="skip cutouts", variable=self.cd_nocut).grid(row=5, column=0, columnspan=2, sticky="w", padx=4)
        self.run_buttons(f, self.run_cand)

    def pick_cand_dir(self):
        d = filedialog.askdirectory(initialdir=self.cd_out.get())
        if d:
            self.cd_out.set(d)

    def run_cand(self):
        try:
            trange = self.cd_int.trange()
            ecut = str(float(self.cd_ecut.get()))
        except ValueError as e:
            messagebox.showerror("bad input", str(e))
            return
        probe, prof, out = self.cd_int.probe.get(), self.cd_profile.get(), Path(self.cd_out.get().strip())
        cmd = [sys.executable, "-u", SCRIPTS / "generate_candidates.py", "--probe", probe, "--trange", *trange,
               "--tag", self.cd_tag.get(), "--profile", prof, "--energy-cutoff", ecut,
               "--out", out, "--datatype", self.cd_dtype.get()]
        if self.cd_nocut.get():
            cmd.append("--no-cutouts")
        jsonl = out / f"{probe}_{trange[0].replace('/', 'T').replace(':', '')}_{prof}.jsonl"
        self.runner.start(cmd, self.env(), lambda rc: self.results.load_candidates(jsonl))

    def close(self):
        self.runner.cancel()
        self.destroy()


if __name__ == "__main__":
    App().mainloop()
