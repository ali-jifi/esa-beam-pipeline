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
import tkinter.font as tkfont
import webbrowser
from datetime import date, datetime, timedelta
from pathlib import Path
from tkinter import filedialog, messagebox, simpledialog, ttk

from PIL import Image, ImageDraw, ImageTk

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


# windows xp luna palette
FACE, WHITE, INK, SEL, TIP = "#ECE9D8", "#FFFFFF", "#000000", "#316AC5", "#FFFFE1"
PANE, PANEL, HEAD_INK, FIELD_LINE = "#6F8CDE", "#D6DFF7", "#215DC6", "#7F9DB9"
F = {}


def _grad(w, h, a, b, horiz=False):
    m = Image.linear_gradient("L")
    m = (m.rotate(90) if horiz else m).resize((w, h))
    return Image.composite(Image.new("RGB", (w, h), b), Image.new("RGB", (w, h), a), m)


def _shape(w, h, r, top, bot, line, outline, inset=None, open_bottom=False, horiz=False):
    # rounded gradient box, open_bottom drops the bottom edge for tabs
    hh = h + r + 2 * line if open_bottom else h
    img = Image.new("RGBA", (w, hh), (0, 0, 0, 0))
    mask = Image.new("L", (w, hh), 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, w - 1, hh - 1), r, fill=255)
    img.paste(_grad(w, hh, top, bot, horiz), (0, 0), mask)
    d = ImageDraw.Draw(img)
    d.rounded_rectangle((0, 0, w - 1, hh - 1), r, outline=outline, width=line)
    if inset:
        d.rounded_rectangle((line, line, w - 1 - line, hh - 1 - line), max(r - line, 0), outline=inset, width=2 * line)
    return img.crop((0, 0, w, h))


def panel_header(w, h, s):
    # task pane header, white to pale blue, rounded top corners on the pane
    r = round(3 * s)
    img = Image.new("RGB", (w, h), PANE)
    mask = Image.new("L", (w, h), 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, w - 1, h + r), r, fill=255)
    img.paste(_grad(w, h, "#FFFFFF", "#C6D3F7", horiz=True), (0, 0), mask)
    return img


def apply_theme(root):
    s = root.winfo_fpixels("1i") / 96
    px = lambda n: max(1, round(n * s))
    line = max(1, int(s))
    F.update(ui="Tahoma", mono="Lucida Console", scale=s)
    root.configure(bg=FACE)
    root.option_add("*Font", (F["ui"], 8))
    root.option_add("*Text.font", (F["mono"], 8))
    root.option_add("*TCombobox*Listbox.font", (F["ui"], 8))
    root.option_add("*TCombobox*Listbox.background", WHITE)
    root.option_add("*TCombobox*Listbox.selectBackground", SEL)
    root.option_add("*TCombobox*Listbox.selectForeground", WHITE)
    st = ttk.Style(root)
    st.theme_use("clam")
    st.configure(".", background=FACE, foreground=INK, bordercolor=FIELD_LINE, lightcolor=WHITE, darkcolor=WHITE,
                 troughcolor=FACE, fieldbackground=WHITE, selectbackground=SEL, selectforeground=WHITE,
                 font=(F["ui"], 8))
    imgs = root._xp_imgs = {}

    def el(name, layout_style, images, layout, **kw):
        ph = [ImageTk.PhotoImage(i, master=root) for i in images[1]]
        imgs[name] = ph
        st.element_create(name, "image", ph[0], *[(spec, p) for spec, p in zip(images[0], ph[1:])], **kw)
        if layout_style:
            st.layout(layout_style, layout)

    w, r = px(21), px(3)
    el("XPButton.border", "TButton",
       (["disabled", "pressed", "active", "focus"],
        [_shape(w, w, r, "#FFFFFF", "#E6E3D3", line, "#003C74"),
         _shape(w, w, r, "#F5F4EA", "#F5F4EA", line, "#C9C7BA"),
         _shape(w, w, r, "#E2E1D9", "#F7F6F2", line, "#003C74"),
         _shape(w, w, r, "#FFFFFF", "#E6E3D3", line, "#003C74", inset="#FAC55E"),
         _shape(w, w, r, "#FFFFFF", "#E6E3D3", line, "#003C74", inset="#98B4EB")]),
       [("XPButton.border", {"sticky": "nswe", "children": [
           ("Button.padding", {"sticky": "nswe", "children": [("Button.label", {"sticky": "nswe"})]})]})],
       border=px(5), sticky="nswe")
    st.configure("TButton", padding=(px(8), px(2)), anchor="center")
    st.map("TButton", foreground=[("disabled", "#A19F94")])

    tw, th = px(20), px(18)

    def tab(fill, orange):
        im = _shape(tw, th, px(3), WHITE, fill, line, "#919B9C", open_bottom=True)
        if orange:
            d = ImageDraw.Draw(im)
            d.rectangle((px(2), 0, tw - 1 - px(2), line - 1), fill="#E68B2C")
            d.rectangle((line, line, tw - 1 - line, 3 * line - 1), fill="#FFC73C")
        return im
    el("XPTab.tab", "TNotebook.Tab",
       (["selected", "active"], [tab("#ECEBE6", False), tab("#FCFCFE", True), tab("#ECEBE6", True)]),
       [("XPTab.tab", {"sticky": "nswe", "children": [("Notebook.padding", {"side": "top", "sticky": "nswe", "children": [
           ("Notebook.label", {"side": "top", "sticky": ""})]})]})],
       border=(px(4), px(4), px(4), 0), sticky="nswe")
    st.configure("TNotebook", bordercolor="#919B9C", lightcolor="#FCFCFE", darkcolor="#FCFCFE", tabmargins=(px(2), px(2), px(2), 0))
    st.configure("TNotebook.Tab", padding=(px(10), px(3)))
    st.map("TNotebook.Tab", expand=[("selected", (px(2), px(2), px(2), 0))])

    b = px(13)

    def box(hover, check):
        im = Image.new("RGBA", (b + px(6), b), (0, 0, 0, 0))
        im.paste(_grad(b, b, "#DCDCD7", WHITE), (0, 0))
        d = ImageDraw.Draw(im)
        d.rectangle((0, 0, b - 1, b - 1), outline="#1C5180", width=line)
        if hover:
            d.rectangle((line, line, b - 1 - line, b - 1 - line), outline="#FAC55E", width=2 * line)
        if check:
            d.line([(px(3), px(6)), (px(5), px(8)), (px(9), px(4))], fill="#21A121", width=2 * line)
        return im
    el("XPCheck.indicator", None,
       (["selected active", "selected", "active"], [box(False, False), box(True, True), box(False, True), box(True, False)]),
       None, sticky="")
    st.layout("TCheckbutton", [("Checkbutton.padding", {"sticky": "nswe", "children": [
        ("XPCheck.indicator", {"side": "left", "sticky": ""}),
        ("Checkbutton.focus", {"side": "left", "sticky": "w", "children": [("Checkbutton.label", {"sticky": "nswe"})]})]})])
    for cls in ("TCheckbutton", "In.TCheckbutton", "Pane.TCheckbutton"):
        st.layout(cls, st.layout("TCheckbutton"))

    ah = tkfont.Font(root, font=(F["ui"], 8)).metrics("linespace") + px(2)
    aw = px(17)

    def arrow(top, bot):
        im = _shape(aw, ah, px(2), top, bot, line, "#B4C8F6")
        cx, cy = aw // 2, ah // 2
        ImageDraw.Draw(im).line([(cx - px(3), cy - px(1)), (cx, cy + px(2)), (cx + px(3), cy - px(1))],
                                fill="#4D6185", width=2 * line)
        return im
    el("XPCombo.arrow", None,
       (["pressed", "active"], [arrow("#E1EAFE", "#B3C7F5"), arrow("#B8CBF6", "#E1EAFE"), arrow("#FDFFFF", "#C2D6FD")]),
       None, border=px(3), sticky="ns")
    st.layout("TCombobox", [("Combobox.field", {"sticky": "nswe", "children": [
        ("XPCombo.arrow", {"side": "right", "sticky": "ns"}),
        ("Combobox.padding", {"sticky": "nswe", "children": [("Combobox.textarea", {"sticky": "nswe"})]})]})])
    st.configure("TCombobox", padding=(px(2), 0))
    st.map("TCombobox", fieldbackground=[("readonly", WHITE)], selectbackground=[("readonly", WHITE)],
           selectforeground=[("readonly", INK)])
    st.configure("TEntry", padding=(px(2), 0), insertcolor=INK)

    hw, hh = px(12), px(18)

    def head(fill, lines):
        im = Image.new("RGB", (hw, hh), fill)
        d = ImageDraw.Draw(im)
        for i, c in enumerate(lines):
            d.rectangle((0, hh - (3 - i) * line, hw - 1, hh - (2 - i) * line - 1), fill=c)
        d.line([(hw - 2 * line, px(3)), (hw - 2 * line, hh - px(5))], fill="#C7C5B2", width=line)
        d.line([(hw - line, px(3)), (hw - line, hh - px(5))], fill=WHITE, width=line)
        return im
    el("XPHead.cell", None,
       (["active"], [head("#EBEADB", ["#E2DECD", "#D6D2C2", "#CBC7B8"]),
                     head("#FAF9F4", ["#F9B119", "#F8A900", "#E59700"])]),
       None, border=(px(2), px(2), px(3), px(4)), sticky="nswe")
    st.layout("Treeview.Heading", [("XPHead.cell", {"sticky": "nswe", "children": [
        ("Treeheading.padding", {"sticky": "nswe", "children": [("Treeheading.text", {"sticky": "w"})]})]})])
    st.configure("Treeview.Heading", padding=(px(5), px(2)))
    st.configure("Treeview", background=WHITE, fieldbackground=WHITE, rowheight=px(17), bordercolor=FIELD_LINE)
    st.map("Treeview", background=[("selected", SEL)], foreground=[("selected", WHITE)])

    for o in ("Vertical", "Horizontal"):
        st.configure(f"{o}.TScrollbar", background="#C1D3FB", troughcolor="#F7F6F0", bordercolor="#98B1E6",
                     lightcolor="#E1EAFE", darkcolor="#A8BEF5", arrowcolor="#4D6185", gripcount=4, arrowsize=px(15))
        st.map(f"{o}.TScrollbar", background=[("pressed", "#A8BEF5"), ("active", "#D7E3FE")])
    st.configure("Sash", sashthickness=px(5), gripcount=0)

    st.configure("Pane.TFrame", background=PANE)
    st.configure("Pane.TLabel", background=PANE, foreground=WHITE, font=(F["ui"], 8, "bold"))
    st.configure("Pane.TCheckbutton", background=PANE, foreground=WHITE, font=(F["ui"], 8, "bold"))
    st.configure("Panel.TFrame", background=PANEL, relief="solid", borderwidth=line, bordercolor=WHITE)
    for cls in ("TFrame", "TLabel", "TCheckbutton"):
        st.configure(f"In.{cls}", background=PANEL)
    for cls in ("Pane", "In"):
        bg = PANE if cls == "Pane" else PANEL
        st.map(f"{cls}.TCheckbutton", background=[("active", bg)])


def walk(w):
    for c in w.winfo_children():
        yield c
        yield from walk(c)


class Panel(ttk.Frame):
    # xp task pane box, gradient header with bold blue title, children grid into the body
    def __init__(self, master, text, **kw):
        font = tkfont.Font(master, font=(F["ui"], 8, "bold"))
        self.head_h = font.metrics("linespace") + round(10 * F["scale"])
        pad = round(8 * F["scale"])
        top = self.head_h + pad // 2
        super().__init__(master, style="Panel.TFrame", padding=(pad, top, pad, pad // 2), **kw)
        self.caption, self._font, self._img = text[:1].upper() + text[1:], font, None
        self.head = tk.Canvas(self, height=self.head_h, highlightthickness=0, bd=0, bg=PANE)
        # place is measured from inside the ttk padding, so step back out to the frame edge
        self.head.place(x=-pad, y=-top, relwidth=1, width=2 * pad)
        self.head.bind("<Configure>", self._draw)

    def _draw(self, e):
        self._img = ImageTk.PhotoImage(panel_header(max(e.width, 1), self.head_h, F["scale"]), master=self)
        self.head.delete("all")
        self.head.create_image(0, 0, image=self._img, anchor="nw")
        self.head.create_text(round(12 * F["scale"]), self.head_h // 2, text=self.caption, anchor="w",
                              font=self._font, fill=HEAD_INK)


def skin(root):
    # restyle plain frames, labels, and checkboxes by where they sit, pane or panel body
    cap = lambda s: s[:1].upper() + s[1:]
    for w in walk(root):
        if not isinstance(w, ttk.Widget):
            continue
        cls = w.winfo_class()
        if cls in ("TButton", "TCheckbutton") and len(str(w.cget("text"))) > 1:
            w.configure(text=cap(str(w.cget("text"))))
        if cls == "TLabel" and len(str(w.cget("text"))) > 2:
            w.configure(text=cap(str(w.cget("text"))))
        if cls not in ("TFrame", "TLabel", "TCheckbutton") or str(w.cget("style")):
            continue
        p, ctx = w.master, None
        while p is not None and ctx is None:
            if isinstance(p, Panel):
                ctx = "In"
            elif isinstance(p, ttk.Frame) and str(p.cget("style")) == "Pane.TFrame":
                ctx = "Pane"
            p = p.master
        if ctx:
            w.configure(style=f"{ctx}.{cls}")


def text_box(parent, **kw):
    # xp edit box, thin blue-gray rule, themed scrollbar
    f = tk.Frame(parent, bg=WHITE, highlightthickness=1, highlightbackground=FIELD_LINE, highlightcolor=FIELD_LINE)
    t = tk.Text(f, relief="flat", bd=0, bg=WHITE, highlightthickness=0, padx=4, pady=2, **kw)
    sb = ttk.Scrollbar(f, orient="vertical", command=t.yview)
    t.configure(yscrollcommand=sb.set)
    sb.pack(side="right", fill="y")
    t.pack(side="left", fill="both", expand=True)
    return f, t


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
        tk.Label(self.top, text=self.text, bg=TIP, fg=INK, font=(F["ui"], 8), relief="solid", bd=1,
                 padx=3, pady=1).pack()
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
        self.configure(bg=FACE)
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
    ttk.Label(parent, text=label[:1].upper() + label[1:] + ":").grid(row=r, column=2 * col, sticky="w", padx=4, pady=1)
    widget.grid(row=r, column=2 * col + 1, sticky="w", padx=4, pady=1)
    if tip:
        Tip(widget, tip)


def date_field(parent, var):
    f = ttk.Frame(parent)
    ttk.Entry(f, textvariable=var, width=11).pack(side="left")
    ttk.Button(f, text="...", width=3, command=lambda: CalendarPopup(parent, var)).pack(side="left")
    return f


class IntervalPicker(Panel):
    # probe + start/end date and time, or a known interval from an events csv
    def __init__(self, master, events: dict):
        super().__init__(master, "interval")
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
    def __init__(self, root, log: tk.Text, status: tk.StringVar):
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
        super().__init__(master, bg=WHITE, highlightthickness=1, highlightbackground=FIELD_LINE)
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
        # buttons pack first so a narrow window shrinks the dropdowns, not the buttons
        ttk.Label(top, text="previous runs", width=12).pack(side="left")
        ttk.Button(top, text="refresh", command=self.refresh_runs).pack(side="right")
        self.runs = ttk.Combobox(top, state="readonly", width=10)
        self.runs.pack(side="left", padx=4, fill="x", expand=True)
        self.runs.bind("<<ComboboxSelected>>", lambda e: self.load_prefix(self.run_paths[self.runs.current()]))
        mid = ttk.Frame(self)
        mid.pack(fill="x", pady=2)
        ttk.Label(mid, text="figure", width=12).pack(side="left")
        ttk.Button(mid, text="fit", command=lambda: self.view.fit()).pack(side="right")
        ttk.Button(mid, text="open folder", command=lambda: self._open(True)).pack(side="right", padx=4)
        ttk.Button(mid, text="open", command=lambda: self._open(False)).pack(side="right")
        self.fig = ttk.Combobox(mid, state="readonly", width=10)
        self.fig.pack(side="left", padx=4, fill="x", expand=True)
        self.fig.bind("<<ComboboxSelected>>", self.show_fig)
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
        tabs.add(tf, text="Table")
        self.tree = ttk.Treeview(tf, show="headings")
        ys = ttk.Scrollbar(tf, orient="vertical", command=self.tree.yview)
        xs = ttk.Scrollbar(tf, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        ys.pack(side="right", fill="y")
        xs.pack(side="bottom", fill="x")
        self.tree.pack(fill="both", expand=True)
        df, self.text = text_box(tabs, height=8)
        tabs.add(df, text="Diagnose")

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
        font = tkfont.Font(self, font=ttk.Style().lookup("Treeview.Heading", "font") or (F["ui"], 8))
        for c in cols:
            self.tree.heading(c, text=c)
            # text plus heading padding and cell border
            self.tree.column(c, width=max(70, font.measure(c) + round(16 * F["scale"])), stretch=False)
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
        apply_theme(self)
        self.defaults = classifier_defaults()
        self.events = load_events()
        self.presets = builtin_presets(self.defaults) | load_user_presets()
        self.data_dir = tk.StringVar(value=os.environ.get("THM_DATA_DIR") or str(ROOT / "data"))
        self.status = tk.StringVar(value="idle")
        self._page_syncs, self._scrollers = [], []

        # xp address bar
        top = ttk.Frame(self, padding=(6, 4))
        top.pack(fill="x")
        ttk.Label(top, text="Data dir", foreground="#7F7C70").pack(side="left")
        e = ttk.Entry(top, textvariable=self.data_dir, width=80)
        e.pack(side="left", padx=4, fill="x", expand=True)
        Tip(e, "where themis cdfs live, missing files download from the themis server on first run")
        ttk.Button(top, text="browse", command=self.pick_data_dir).pack(side="left")
        # xp status bar
        bar = tk.Frame(self, bg=FACE)
        bar.pack(side="bottom", fill="x")
        tk.Label(bar, textvariable=self.status, bg=FACE, anchor="w", relief="sunken", bd=1, padx=4,
                 width=30).pack(side="left", padx=(2, 0), pady=2)
        tk.Label(bar, text="Figure: wheel zooms, drag pans, double click refits", bg=FACE, anchor="w",
                 relief="sunken", bd=1, padx=4).pack(
            side="left", fill="x", expand=True, padx=2, pady=2)

        vpane = ttk.PanedWindow(self, orient="vertical")
        vpane.pack(fill="both", expand=True)
        vpane.bind("<Configure>", lambda e: (vpane.sashpos(0, int(e.height * 0.8)), vpane.unbind("<Configure>"))
                   if e.height > 100 else None)
        hpane = ttk.PanedWindow(vpane, orient="horizontal")
        vpane.add(hpane, weight=4)
        self.tabs = ttk.Notebook(hpane)
        hpane.add(self.tabs, weight=0)
        self.results = ResultsPanel(hpane)
        hpane.add(self.results, weight=1)
        lf, self.log = text_box(vpane, height=10)
        vpane.add(lf, weight=1)
        self.runner = Runner(self, self.log, self.status)

        self.build_pipeline_tab()
        self.build_spectro_tab()
        self.build_cand_tab()
        skin(self)
        self.update_idletasks()
        for sync in self._page_syncs:
            sync()
        self.minsize(960, 600)
        self.bind_all("<MouseWheel>", self._wheel, add="+")
        for t in self.tabs.tabs():
            txt = self.tabs.tab(t, "text")
            self.tabs.tab(t, text=txt[:1].upper() + txt[1:])
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
        bf.pack(fill="x", pady=10)
        ttk.Button(bf, text="run", command=cmd).pack(side="left")
        ttk.Button(bf, text="cancel", command=self.runner.cancel).pack(side="left", padx=4)

    def build_pipeline_tab(self):
        f = self.page("beam pipeline")
        self.pp_int = IntervalPicker(f, self.events)
        self.pp_int.pack(fill="x")

        pf = ttk.Frame(f)
        pf.pack(fill="x", pady=8)
        ttk.Label(pf, text="preset").pack(side="left")
        self.preset = ttk.Combobox(pf, values=list(self.presets), state="readonly", width=18)
        self.preset.pack(side="left", padx=4)
        self.preset.bind("<<ComboboxSelected>>", lambda e: self.apply_preset(self.preset.get()))
        ttk.Button(pf, text="save as", command=self.save_preset).pack(side="left")

        self.pvars = {k: tk.StringVar() for k, *_ in EXTRA + PARAMS}
        basic = Panel(f, "parameters")
        basic.pack(fill="x")
        for r, (key, _, label, tip, _) in enumerate(EXTRA + [p for p in PARAMS if p[4]]):
            row(basic, r, label, ttk.Entry(basic, textvariable=self.pvars[key], width=8), tip)

        self.show_adv = tk.BooleanVar()
        ttk.Checkbutton(f, text="show advanced", variable=self.show_adv, command=self.toggle_adv).pack(anchor="w", pady=2)
        self.adv = Panel(f, "advanced")
        for i, (key, _, label, tip, _) in enumerate([p for p in PARAMS if not p[4]]):
            row(self.adv, i // 2, label, ttk.Entry(self.adv, textvariable=self.pvars[key], width=8), tip, col=i % 2)

        self.opt = Panel(f, "options")
        self.opt.pack(fill="x", pady=(8, 0))
        self.tcv = tk.StringVar(value="1.0 1.2 1.5 2.0")
        row(self.opt, 0, "threshold compare R", ttk.Entry(self.opt, textvariable=self.tcv, width=20),
            "coherent R values for the threshold-comparison plot, space separated")
        self.no_plots = tk.BooleanVar()
        ttk.Checkbutton(self.opt, text="skip plots", variable=self.no_plots).grid(row=1, column=0, sticky="w", padx=4)
        self.diag = tk.BooleanVar()
        ttk.Checkbutton(self.opt, text="diagnose dump", variable=self.diag).grid(row=1, column=1, sticky="w", padx=4)
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
            self.adv.pack(fill="x", pady=(8, 0), before=self.opt)
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
        f = self.page("spectrogram")
        g = Panel(f, "ion eflux")
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
        f = self.page("candidates")
        self.cd_int = IntervalPicker(f, self.events)
        self.cd_int.pack(fill="x")
        g = Panel(f, "generation")
        g.pack(fill="x", pady=(10, 0))
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

    def page(self, text):
        # task pane page that scrolls when the window is too short for it
        outer = ttk.Frame(self.tabs, style="Pane.TFrame")
        self.tabs.add(outer, text=text)
        cv = tk.Canvas(outer, bg=PANE, highlightthickness=0, bd=0, height=200, yscrollincrement=20)
        sb = ttk.Scrollbar(outer, orient="vertical", command=cv.yview)
        cv.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        cv.pack(side="left", fill="both", expand=True)
        inner = ttk.Frame(cv, padding=10, style="Pane.TFrame")
        win = cv.create_window(0, 0, window=inner, anchor="nw")

        def sync(_=None):
            cv.configure(width=inner.winfo_reqwidth(), scrollregion=(0, 0, inner.winfo_reqwidth(), inner.winfo_reqheight()))
        inner.bind("<Configure>", sync)
        cv.bind("<Configure>", lambda e: (cv.itemconfigure(win, width=e.width), sync()))
        self._page_syncs.append(sync)
        self._scrollers.append(cv)
        return inner

    def _wheel(self, e):
        # wheel over a task pane scrolls it, combos keep their own wheel handling
        try:
            w = e.widget if not isinstance(e.widget, str) else self.nametowidget(e.widget)
        except (KeyError, tk.TclError):
            return
        if isinstance(w, ttk.Combobox):
            return
        while w is not None:
            if w in self._scrollers:
                top, bot = w.yview()
                if top > 0 or bot < 1:
                    w.yview_scroll(-1 if e.delta > 0 else 1, "units")
                return
            w = w.master

    def close(self):
        self.runner.cancel()
        self.destroy()


if __name__ == "__main__":
    App().mainloop()
