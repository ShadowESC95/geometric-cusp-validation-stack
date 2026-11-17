"""Desktop window for the cusp functional solver.

A Tk window with embedded Matplotlib figures. All numbers come from
``geocusp.solver`` in this process. There is no server and no network access.
"""
from __future__ import annotations

import tkinter as tk
from math import log10, pi
from tkinter import font as tkfont
from tkinter import ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from . import solver

BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#898781"
GOOD, BAD = "#006300", "#d03b3b"
MONO = "TkFixedFont"  # Tk's own fixed-width font, present on every platform
REGIMES = {
    "collision": "Double root on the rim. A single-root residue grows like |disc|^(-5/2) here and "
                 "its absolute integral diverges logarithmically; the paired period stays finite.",
    "one-meridian": "One root inside the rim. This is the one-meridian band, where a single "
                    "residue survives at first order in the regulator.",
    "paired": "Both roots inside the rim. The rim contour encloses both meridians, "
              "so the paired period is what it measures.",
    "outside": "Both roots outside the rim. Neither meridian contributes on this slice.",
}


def signed(value: float, digits: int = 9) -> str:
    return f"{value:+.{digits}f}"


def cplx(pair, digits: int = 4) -> str:
    return f"{pair[0]:+.{digits}f} {'-' if pair[1] < 0 else '+'} {abs(pair[1]):.{digits}f}i"


def cplx_e(pair) -> str:
    return f"{pair[0]:+.3e} {'-' if pair[1] < 0 else '+'} {abs(pair[1]):.3e}i"


class Control:
    """A labelled slider with a typed entry kept in sync."""

    def __init__(self, parent, row, label, low, high, value, on_change, digits=4,
                 integer=False, log=False):
        self.low, self.high, self.digits = low, high, digits
        self.integer, self.log, self.on_change = integer, log, on_change
        self.exact = None
        self.text = tk.StringVar()
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
        self.scale = ttk.Scale(parent, from_=self._to_slider(low), to=self._to_slider(high),
                               command=self._moved)
        self.scale.grid(row=row, column=1, sticky="ew", pady=2)
        self.entry = ttk.Entry(parent, textvariable=self.text, width=10, justify="right", font=MONO)
        self.entry.grid(row=row, column=2, padx=(8, 0), pady=2)
        self.entry.bind("<Return>", self._typed)
        self.entry.bind("<FocusOut>", self._typed)
        self.set(value, notify=False)

    def _to_slider(self, value):
        return log10(value) if self.log else value

    def _show(self, value):
        self.text.set(str(int(value)) if self.integer else f"{value:.{self.digits}f}")

    def _moved(self, _):
        raw = self.scale.get()
        value = 10**raw if self.log else raw
        self.exact = round(value) if self.integer else value
        self._show(self.exact)
        self.on_change()

    def _typed(self, _):
        try:
            value = float(self.text.get())
        except ValueError:
            self._show(self.get())
            return
        if value != self.get():
            self.set(value)

    def set(self, value, notify=True):
        value = min(self.high, max(self.low, value))
        self.exact = round(value) if self.integer else value
        self.scale.configure(command="")
        self.scale.set(self._to_slider(self.exact))
        self.scale.configure(command=self._moved)
        self._show(self.exact)
        if notify:
            self.on_change()

    def get(self):
        return self.exact


class Table(ttk.Treeview):
    def __init__(self, parent, columns, widths, height):
        super().__init__(parent, columns=columns, show="headings", height=height)
        for name, width in zip(columns, widths):
            self.heading(name, text=name, anchor="w")
            self.column(name, width=width, anchor="w", stretch=True)
        self.tag_configure("ok", foreground=GOOD)
        self.tag_configure("bad", foreground=BAD)

    def fill(self, rows):
        self.delete(*self.get_children())
        for *values, tag in rows:
            self.insert("", "end", values=values, tags=(tag,))


class SolverApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("Cusp Functional Solver")
        width = min(1280, root.winfo_screenwidth()-40)
        height = min(780, root.winfo_screenheight()-80)
        root.geometry(f"{width}x{height}")
        root.minsize(min(900, width), min(620, height))
        self.pending = {}
        self.results = {}
        style = ttk.Style(root)
        self.big_font = tkfont.nametofont(MONO).copy()
        self.big_font.configure(size=16)
        style.configure("Big.TLabel", font=self.big_font)
        style.configure("Cap.TLabel", foreground=GREY, font=("TkDefaultFont", 8))
        style.configure("Treeview", font=MONO, rowheight=24)
        book = ttk.Notebook(root)
        book.pack(fill="both", expand=True, padx=8, pady=8)
        for title, build in (("Functional solver", self.build_functional),
                             ("Lorentzian edge chain", self.build_edge),
                             ("Integration contour", self.build_contour),
                             ("CD weight at a point", self.build_cd),
                             ("Rim collision explorer", self.build_rim)):
            page = ttk.Frame(book, padding=12)
            book.add(page, text=title)
            page.columnconfigure(1, weight=1)
            page.rowconfigure(0, weight=1)
            build(page)
        self.update_functional()
        self.update_edge()
        self.update_cd()
        self.update_rim()

    # -- shared helpers

    def later(self, key, action, delay=120):
        if key in self.pending:
            self.root.after_cancel(self.pending[key])
        self.pending[key] = self.root.after(delay, lambda: (self.pending.pop(key, None), action()))

    def panel(self, page):
        left = ttk.Frame(page)
        left.grid(row=0, column=0, sticky="nsw", padx=(0, 16))
        left.columnconfigure(1, weight=1, minsize=190)
        right = ttk.Frame(page)
        right.grid(row=0, column=1, sticky="nsew")
        right.columnconfigure(0, weight=1)
        return left, right

    def figure(self, parent, row, size):
        figure = Figure(figsize=size, dpi=100, layout="constrained")
        canvas = FigureCanvasTkAgg(figure, master=parent)
        canvas.get_tk_widget().grid(row=row, column=0, sticky="nsew", pady=(8, 0))
        parent.rowconfigure(row, weight=1)
        return figure, canvas

    def note(self, parent, row):
        label = ttk.Label(parent, wraplength=640, justify="left")
        label.grid(row=row, column=0, sticky="w", pady=(8, 0))
        return label

    # -- functional solver

    def build_functional(self, page):
        left, right = self.panel(page)
        self.formula = ttk.Label(left, font=MONO, wraplength=380, justify="left")
        self.formula.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
        buttons = ttk.Frame(left)
        buttons.grid(row=1, column=0, columnspan=3, sticky="w", pady=(0, 8))
        presets = (("F0 = -1", (-1, 0, 0, 1)), ("F1 = log²z + π²", (pi**2, 1, 0, 2)),
                   ("z^p only", (0, 0, 1, 2)))
        for column, (text, values) in enumerate(presets):
            ttk.Button(buttons, text=text, command=lambda v=values: self.preset(*v)).grid(
                row=0, column=column, padx=(0, 6))
        change = lambda: self.later("functional", self.update_functional)
        specs = (("c0", "constant c0", -12, 12, pi**2), ("c1", "c1 · log z", -3, 3, 0),
                 ("c2", "c2 · log² z", -3, 3, 1), ("c3", "c3 · log³ z", -3, 3, 0),
                 ("c4", "c4 · log⁴ z", -3, 3, 0), ("A", "A · z^p", -3, 3, 0),
                 ("p", "exponent p", -.95, .95, .3))
        self.f = {key: Control(left, row, label, low, high, value, change)
                  for row, (key, label, low, high, value) in enumerate(specs, start=2)}
        self.f["L"] = Control(left, 9, "loop order L", 1, 4, 2, change, integer=True)
        self.f["R"] = Control(left, 10, "cycle radius |a|", .3, 3, 1.3, change, digits=2)
        self.f["eta"] = Control(left, 11, "regulator -eps", .0025, .1, .01, change, log=True)

        heads = ttk.Frame(right)
        heads.grid(row=0, column=0, sticky="w")
        self.big = {}
        for column, (key, caption) in enumerate((("mean", "MEAN M[F]"), ("gamma", "GAMMA(L) = -4M/L"),
                                                 ("reference", "REFERENCE"))):
            ttk.Label(heads, text=caption, style="Cap.TLabel").grid(row=0, column=column, sticky="w", padx=(0, 28))
            self.big[key] = ttk.Label(heads, text="-", style="Big.TLabel")
            self.big[key].grid(row=1, column=column, sticky="w", padx=(0, 28))
        self.f_table = Table(right, ("Functional", "How it is computed", "Value", "- closed form"),
                             (190, 330, 140, 110), 4)
        self.f_table.grid(row=1, column=0, sticky="ew", pady=(10, 0))
        self.f_figure, self.f_canvas = self.figure(right, 2, (6.4, 2.8))
        self.verdict = self.note(right, 3)

    def preset(self, c0, c2, amplitude, loop_order):
        for key, value in (("c0", c0), ("c1", 0), ("c2", c2), ("c3", 0), ("c4", 0),
                           ("A", amplitude), ("L", loop_order)):
            self.f[key].set(value, notify=False)
        self.update_functional()

    def weight(self):
        return solver.Weight(tuple(self.f[f"c{k}"].get() for k in range(5)),
                             self.f["A"].get(), self.f["p"].get())

    def update_functional(self):
        weight = self.weight()
        self.formula.configure(text="F(z) = "+weight.describe())
        try:
            r = solver.solve_functional(weight, self.f["L"].get(), self.f["R"].get(),
                                        self.f["eta"].get(), curve=True)
        except (ValueError, ArithmeticError) as error:
            self.verdict.configure(text=f"Solver error: {error}", foreground=BAD)
            return
        self.results["functional"] = r
        exact, scale = r["closed_mean"], max(1., abs(r["closed_mean"]))
        self.big["mean"].configure(text=signed(exact, 6))
        self.big["gamma"].configure(text=signed(r["gamma"], 6))
        known = r["reference"] is not None
        self.big["reference"].configure(
            text=signed(r["reference"], 6) if known else "none",
            foreground=GOOD if r["matches_reference"] and r["loop_order"] <= 2 else "")
        rows = (("Circle mean", "Gauss-Legendre on |z| = 1", r["circle_mean"], 1e-9),
                ("Capped-cycle period", f"3-form on the cycle, |a| = {r['radius']:.2f}",
                 r["cycle_period"], 1e-5),
                ("Lorentzian pole", f"Mellin-Barnes at eps = -{r['eta']:.4f}", r["lorentzian"], None),
                ("Lorentzian pole, eps -> 0", "Richardson step from eps and eps/2",
                 r["lorentzian_extrapolated"], 2e-3))
        self.f_table.fill([(name, how, signed(value, 8), f"{value-exact:+.1e}",
                            "" if tol is None else ("ok" if abs(value-exact) < tol*scale else "bad"))
                           for name, how, value, tol in rows])
        text = ("All three functionals agree on this weight. " if r["agree"] else
                "The functionals disagree beyond tolerance; lower the regulator or check the inputs. ")
        if not known:
            text += "No reference coefficient is built in at this loop order."
        elif r["loop_order"] == 3:
            text += ("The three-loop coefficient comes from the two-loop weight F2, which is not "
                     "in this polynomial family. Use the Lorentzian edge chain tab.")
        elif r["matches_reference"]:
            text += f"With L = {r['loop_order']} this weight reproduces the known coefficient."
        else:
            text += f"This weight does not give the known L = {r['loop_order']} coefficient."
        self.verdict.configure(text=text, foreground=GOOD if r["agree"] else BAD)

        self.f_figure.clear()
        axes = self.f_figure.add_subplot()
        etas = [point["eta"] for point in r["curve"]]
        axes.plot(etas, [point["value"] for point in r["curve"]], "o-", color=BLUE, lw=2, ms=6,
                  label="eps²·s·t·∫F / 4 at regulator -eps")
        axes.axhline(exact, color=ORANGE, lw=2, ls="--", label="circle mean (closed form)")
        axes.set_xscale("log")
        axes.invert_xaxis()
        axes.set_xlabel("regulator -eps (shrinking to the right)")
        axes.set_title("Lorentzian double pole as the regulator is removed", fontsize=10, loc="left")
        axes.grid(True, color="#e1e0d9", lw=.8)
        axes.legend(fontsize=8, frameon=False)
        self.f_canvas.draw_idle()

    # -- Lorentzian edge chain

    def build_edge(self, page):
        left, right = self.panel(page)
        self.e_formula = ttk.Label(left, font=MONO, wraplength=380, justify="left")
        self.e_formula.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 4))
        ttk.Label(left, text="The weight is the one set on the Functional solver tab.",
                  style="Cap.TLabel", wraplength=380, justify="left").grid(
            row=1, column=0, columnspan=3, sticky="w", pady=(0, 8))
        change = lambda: self.later("edge", self.update_edge, delay=400)
        self.e = {"cusp": Control(left, 2, "cusp", 1, 4, 2, change, integer=True),
                  "position": Control(left, 3, "position on the line", -.9, 1.9, .4, change, digits=3)}
        self.e_edge = tk.StringVar(value="next")
        self.e_region = tk.StringVar(value="--")
        for row, (label, variable, values) in enumerate(
                (("edge at the cusp", self.e_edge, ("next", "prev")),
                 ("signs of s and t", self.e_region, tuple(solver.REGIONS))), start=4):
            ttk.Label(left, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            box = ttk.Combobox(left, textvariable=variable, values=values, state="readonly", width=8)
            box.grid(row=row, column=1, sticky="w", pady=2)
            box.bind("<<ComboboxSelected>>", lambda _: change())
        buttons = ttk.Frame(left)
        buttons.grid(row=6, column=0, columnspan=3, sticky="w", pady=(10, 0))
        ttk.Button(buttons, text="Scan the null line", command=self.scan_edge).grid(
            row=0, column=0, padx=(0, 6))
        ttk.Button(buttons, text="Two-loop weight F2 (about a minute)",
                   command=self.two_loop).grid(row=0, column=1)
        ttk.Label(left, text="The edge segment is 0 < position < 1. Positions 0 and 1 are the "
                             "cusps. On the segment the chain equals the unit-circle mean of "
                             "the weight; on the rest of the null line it is zero.",
                  style="Cap.TLabel", wraplength=380, justify="left").grid(
            row=7, column=0, columnspan=3, sticky="w", pady=(10, 0))

        heads = ttk.Frame(right)
        heads.grid(row=0, column=0, sticky="w")
        self.e_big = {}
        for column, (key, caption) in enumerate((("chain", "EDGE CHAIN  jump / (-2πi)"),
                                                 ("mean", "CIRCLE MEAN M[F]"),
                                                 ("expected", "EXPECTED HERE"))):
            ttk.Label(heads, text=caption, style="Cap.TLabel").grid(row=0, column=column, sticky="w", padx=(0, 28))
            self.e_big[key] = ttk.Label(heads, text="-", style="Big.TLabel")
            self.e_big[key].grid(row=1, column=column, sticky="w", padx=(0, 28))
        self.e_table = Table(right, ("Quantity", "Meaning", "Value"), (200, 400, 170), 5)
        self.e_table.grid(row=1, column=0, sticky="ew", pady=(10, 0))
        self.e_figure, self.e_canvas = self.figure(right, 2, (6.4, 2.8))
        self.e_verdict = self.note(right, 3)

    def edge_settings(self):
        position = self.e["position"].get()
        if abs(position) < 1e-3 or abs(position-1) < 1e-3:
            position += 2e-3          # step off the cusp itself
        return self.e["cusp"].get(), self.e_edge.get(), position, self.e_region.get()

    def update_edge(self):
        weight = self.weight()
        self.e_formula.configure(text="F(z) = "+weight.describe())
        cusp, edge, position, region = self.edge_settings()
        try:
            r = solver.edge_chain(weight, cusp, edge, position, region)
        except (ValueError, ArithmeticError) as error:
            self.e_verdict.configure(text=f"Solver error: {error}", foreground=BAD)
            return
        self.results["edge"] = r
        self.show_edge(r)

    def show_edge(self, r):
        self.e_big["chain"].configure(text=signed(r["edge_chain"], 7))
        self.e_big["mean"].configure(text=signed(r["circle_mean"], 7))
        self.e_big["expected"].configure(text=signed(r["expected"], 7),
                                         foreground=GOOD if r["agree"] else BAD)
        tag = "ok" if r["agree"] else "bad"
        self.e_table.fill([
            ("Edge chain", "jump of the transverse kernel / (-2πi)", signed(r["edge_chain"], 9), tag),
            ("Imaginary part", "should vanish", f"{r['edge_chain_imag']:+.1e}", ""),
            ("Constant weight", "the same with F = 1: 1 on the segment, 0 off it",
             signed(r["constant_weight"], 9), ""),
            ("Sheets visited", "values of arg z / π along the contour",
             ", ".join(str(k) for k in r["sheets"]), ""),
            ("Point", f"cusp {r['cusp']}, {r['edge']} edge, signs {r['region']}",
             "on the segment" if r["on_segment"] else "off the segment", ""),
        ])
        where = "on the edge segment" if r["on_segment"] else "off the edge segment"
        if r["agree"]:
            text = (f"This point is {where}. The Lorentzian chain "
                    + ("reproduces the unit-circle mean." if r["on_segment"]
                       else "vanishes: the contour closes up on one sheet."))
        else:
            text = f"This point is {where}, and the chain does not match the expected value."
        self.e_verdict.configure(text=text, foreground=GOOD if r["agree"] else BAD)

    def scan_edge(self):
        weight = self.weight()
        cusp, edge, _, region = self.edge_settings()
        positions = (-.6, -.3, -.1, .1, .3, .5, .7, .9, 1.1, 1.3, 1.6)
        self.root.configure(cursor="watch")
        self.root.update_idletasks()
        try:
            scan = solver.edge_scan(weight, positions, cusp, edge, region)
        except (ValueError, ArithmeticError) as error:
            self.e_verdict.configure(text=f"Solver error: {error}", foreground=BAD)
            return
        finally:
            self.root.configure(cursor="")
        self.results["edge_scan"] = scan
        self.e_figure.clear()
        axes = self.e_figure.add_subplot()
        axes.axvspan(0, 1, color="#e9f1fb", lw=0)
        axes.plot(positions, [point["edge_chain"] for point in scan], "o", color=BLUE, ms=7,
                  label="edge chain")
        mean = scan[0]["circle_mean"]
        axes.plot([-.7, 0, 0, 1, 1, 1.7], [0, 0, mean, mean, 0, 0], color=ORANGE, lw=2, ls="--",
                  label="circle mean on the segment, zero off it")
        axes.set_xlabel("position along the null line (shaded: the edge segment)")
        axes.set_title("Collinear singularity along the null line through the cusp", fontsize=10, loc="left")
        axes.grid(True, color="#e1e0d9", lw=.8)
        axes.legend(fontsize=8, frameon=False)
        self.e_canvas.draw_idle()
        good = all(point["agree"] for point in scan)
        self.e_verdict.configure(
            text=("All scanned points agree: the mean on the segment, zero off it." if good else
                  "Some scanned points do not match the expected value."),
            foreground=GOOD if good else BAD)

    def two_loop(self):
        self.root.configure(cursor="watch")
        self.e_verdict.configure(text="Computing the two-loop weight on its continued sheets ...",
                                 foreground="")
        self.root.update_idletasks()
        try:
            r = solver.solve_two_loop(edge=True)
        except (ValueError, ArithmeticError) as error:
            self.e_verdict.configure(text=f"Solver error: {error}", foreground=BAD)
            return
        finally:
            self.root.configure(cursor="")
        self.results["two_loop"] = r
        self.e_formula.configure(text="F(z) = two-loop weight F2")
        self.e_big["chain"].configure(text=signed(r["edge_chain"], 7))
        self.e_big["mean"].configure(text=signed(r["circle_mean"], 7))
        self.e_big["expected"].configure(text=signed(r["closed_mean"], 7),
                                         foreground=GOOD if r["agree"] else BAD)
        tag = "ok" if r["agree"] else "bad"
        self.e_table.fill([
            ("Edge chain", "on the segment at cusp 2, continued sheets", signed(r["edge_chain"], 9), tag),
            ("Edge chain off segment", "the same null line beyond the cusp",
             f"{r['edge_chain_off_segment']:+.1e}", ""),
            ("Circle mean", "closed form is -11π⁴/15", signed(r["circle_mean"], 9), ""),
            ("F2 at z = -1", "regular there; closed form 11π⁴/45",
             signed(r["value_at_z_minus_one"], 9), ""),
            ("Γ(3) = -4M/3", f"reference {r['reference']:.9f}", signed(r["gamma"], 9),
             "ok" if r["matches_reference"] else "bad"),
        ])
        self.e_verdict.configure(
            text=("With the two-loop weight the edge chain gives -11π⁴/15, hence the three-loop "
                  "coefficient 44π⁴/45. By Cauchy's theorem this is the same number as the "
                  "circle mean: it tests the continuation of F2 to its second sheet, not an "
                  "independent quantity." if r["agree"] else
                  "The two-loop weight did not reproduce the expected values."),
            foreground=GOOD if r["agree"] else BAD)

    # -- integration contour

    def build_contour(self, page):
        left, right = self.panel(page)
        ttk.Label(left, text="The explicit contour x0 → x0 + iV(x0) of GATE_6_CONTOUR, with "
                             "constants derived from the square.",
                  wraplength=380, justify="left").grid(row=0, column=0, columnspan=3, sticky="w",
                                                       pady=(0, 8))
        self.k_region = tk.StringVar(value="--")
        ttk.Label(left, text="signs of s and t").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Combobox(left, textvariable=self.k_region, values=tuple(solver.REGIONS),
                     state="readonly", width=8).grid(row=1, column=1, sticky="w", pady=2)
        idle = lambda: None
        self.k = {"angle": Control(left, 2, "shape of the square", .5, 2.4, 1.1, idle, digits=2),
                  "points": Control(left, 3, "sample points", 500, 20000, 3000, idle, integer=True),
                  "eta": Control(left, 4, "regulator η", .02, .45, .2, idle, digits=3),
                  "sigma": Control(left, 5, "Re p", -.6, .6, .1, idle, digits=3)}
        buttons = ttk.Frame(left)
        buttons.grid(row=6, column=0, columnspan=3, sticky="w", pady=(10, 0))
        for column, (text, command) in enumerate((("Certificate", self.contour_certificate),
                                                  ("Survey (slow)", self.contour_survey),
                                                  ("Convergence", self.contour_threshold))):
            ttk.Button(buttons, text=text, command=command).grid(row=0, column=column, padx=(0, 6))
        ttk.Label(left, text="Certificate: the inequalities the proof needs, evaluated in floating "
                             "point. Survey: the same properties measured on stress samples, with "
                             "the range of arg z. Convergence: the integral of |z^p| near the edges "
                             "is finite only for |Re p| < η.",
                  style="Cap.TLabel", wraplength=380, justify="left").grid(
            row=7, column=0, columnspan=3, sticky="w", pady=(10, 0))
        self.k_table = Table(right, ("Quantity", "Meaning", "Value"), (230, 400, 170), 18)
        self.k_table.grid(row=0, column=0, sticky="nsew")
        self.k_verdict = self.note(right, 1)
        self.k_verdict.configure(text="Choose a region and press a button.")

    def contour_run(self, message, action):
        self.root.configure(cursor="watch")
        self.k_verdict.configure(text=message, foreground="")
        self.root.update_idletasks()
        try:
            return action()
        except (ValueError, ArithmeticError) as error:
            self.k_verdict.configure(text=f"Solver error: {error}", foreground=BAD)
            return None
        finally:
            self.root.configure(cursor="")

    def contour_certificate(self):
        r = self.contour_run("Evaluating the certificate ...", lambda: solver.contour_certificate(
            self.k_region.get(), self.k["angle"].get()))
        if r is None:
            return
        self.results["certificate"] = r
        rows = [(key, "condition", "holds" if value else "FAILS", "ok" if value else "bad")
                for key, value in r["conditions"].items()]
        rows += [(key, "constant", f"{value:.6g}", "") for key, value in r["constants"].items()]
        self.k_table.fill(rows)
        self.k_verdict.configure(
            text=(f"All {len(r['conditions'])} conditions hold for s = {r['s']:.3f}, t = {r['t']:.3f}. "
                  "The arithmetic is floating point, not interval arithmetic."
                  if r["all"] else "Conditions fail: " + ", ".join(r["failed"])),
            foreground=GOOD if r["all"] else BAD)

    def contour_survey(self):
        r = self.contour_run("Sampling the contour ...", lambda: solver.contour_survey(
            self.k["points"].get(), 0, self.k_region.get(), self.k["angle"].get()))
        if r is None:
            return
        self.results["contour"] = r
        ok = lambda flag: "ok" if flag else "bad"
        low, high = r["arg_z_over_pi_contour"]
        real_low, real_high = r["arg_z_over_pi_real"]
        self.k_table.fill([
            ("Largest Im x0i² near a cone", "must not be positive", f"{r['max_im_near']:+.2e}",
             ok(r["feynman_side_of_every_cone"])),
            ("Near-cone pairs tested", "", str(r["near_count"]), ""),
            ("Smallest |x0i²| / dist²", "zeros only on the edges", f"{r['min_ratio']:.2e}",
             ok(r["zeros_only_on_edges"])),
            ("Excursion of arg x0i²", "beyond [-π, 0], radians", f"{r['max_excursion']:.4f}", ""),
            ("arg z / π at real points", "Lemma P", f"{real_low:+.3f} … {real_high:+.3f}", ""),
            ("arg z / π on the contour", f"proved bound {r['arg_z_over_pi_proved_bound']:.3f}",
             f"{low:+.3f} … {high:+.3f}", ok(r["phase_within_proved_bound"])),
            ("|V²| / |Im x0i²| near cones", f"Lemma Q bound {r['mu_Q']:.4f}",
             f"{r['corner_ratio']:.4f}", ok(r["corner_ratio_within_bound"])),
            ("Certificate", "all conditions", "holds" if r["certificate_passes"] else "FAILS",
             ok(r["certificate_passes"])),
        ])
        good = (r["feynman_side_of_every_cone"] and r["zeros_only_on_edges"]
                and r["phase_within_proved_bound"] and r["corner_ratio_within_bound"])
        self.k_verdict.configure(
            text=("Every sampled point is on the Feynman side of every cone, and arg z stays "
                  "inside the proved bound." if good else
                  "A sampled point violates one of the properties."),
            foreground=GOOD if good else BAD)

    def contour_threshold(self):
        r = self.contour_run("Integrating ...", lambda: solver.convergence_threshold(
            self.k["eta"].get(), self.k["sigma"].get()))
        if r is None:
            return
        self.results["threshold"] = r
        rows = [(f"cutoff {cutoff:.0e}", "strip integral of Lemma 8′", f"{value:.5g}", "")
                for cutoff, value in zip(r["cutoffs"], r["values"])]
        rows.append(("Growth between cutoffs", "settles if convergent",
                     f"{r['growth'][0]:.4g}, {r['growth'][1]:.4g}", ""))
        self.k_table.fill(rows)
        self.k_verdict.configure(
            text=(f"|Re p| = {abs(r['sigma']):.3f} is below η = {r['eta']:.3f}: the integral of "
                  "|z^p| converges, and the exchange of integrals is allowed on this line."
                  if r["converges"] else
                  f"|Re p| = {abs(r['sigma']):.3f} is not below η = {r['eta']:.3f}: the integral "
                  "of |z^p| diverges near the cusps. The proof works on lines inside the "
                  "regulator and moves them afterwards."),
            foreground=GOOD if r["converges"] else BAD)

    # -- CD weight

    def build_cd(self, page):
        left, right = self.panel(page)
        ttk.Label(left, font=MONO, wraplength=380, justify="left",
                  text="F(z) = ∫ 2(1+z)(log y - log z) / ((y - z)(1+y)) dy,  y from 0 to ∞").grid(
            row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
        change = lambda: self.later("cd", self.update_cd, 60)
        self.z = {"re": Control(left, 1, "Re z", -3, 3, .35, change),
                  "im": Control(left, 2, "Im z", -3, 3, .8, change)}
        buttons = ttk.Frame(left)
        buttons.grid(row=3, column=0, columnspan=3, sticky="w", pady=8)
        ttk.Button(buttons, text="Put z on the unit circle", command=self.to_circle).grid(row=0, column=0, padx=(0, 6))
        ttk.Button(buttons, text="Step across the cut", command=self.across_cut).grid(row=0, column=1)
        self.cd_table = Table(right, ("Quantity", "Re", "Im"), (220, 190, 190), 4)
        self.cd_table.grid(row=0, column=0, sticky="ew")
        self.cd_figure, self.cd_canvas = self.figure(right, 1, (4.6, 4.2))
        self.cd_note = self.note(right, 2)

    def to_circle(self):
        x, y = self.z["re"].get(), self.z["im"].get()
        size = (x*x+y*y)**.5 or 1.
        self.z["re"].set(x/size, notify=False)
        self.z["im"].set(y/size)

    def across_cut(self):
        self.z["re"].set(-abs(self.z["re"].get() or 1.), notify=False)
        self.z["im"].set(-.01 if self.z["im"].get() > 0 else .01)

    def update_cd(self):
        z = complex(self.z["re"].get(), self.z["im"].get())
        try:
            r = solver.solve_cd(z)
        except ValueError as error:
            self.cd_table.fill([])
            self.cd_note.configure(text=f"Solver error: {error}", foreground=BAD)
            return
        self.results["cd"] = r
        value, closed = r["value"], r["closed_form"]
        self.cd_table.fill([
            ("Integrated CD weight", signed(value[0]), signed(value[1]), ""),
            ("log²z + π²", signed(closed[0]), signed(closed[1]), ""),
            ("Difference", f"{value[0]-closed[0]:+.1e}", f"{value[1]-closed[1]:+.1e}",
             "ok" if r["error"] < 1e-8 else "bad"),
            ("log z", signed(r["log_z"][0]), signed(r["log_z"][1]), "")])
        if r["on_cut"]:
            text = ("z is on the negative real axis, the principal cut. The value shown is the "
                    "limit from above; from below the imaginary part flips sign.")
        else:
            text = f"|z| = {r['modulus']:.4f}. Quadrature minus closed form: {r['error']:.1e}."
        self.cd_note.configure(text=text, foreground="")
        self.cd_figure.clear()
        axes = self.cd_figure.add_subplot()
        axes.add_patch(_circle(1., GREY))
        axes.plot([-3.2, 0], [0, 0], color=ORANGE, lw=2.5, label="principal cut")
        axes.plot([z.real], [z.imag], "o", color=BLUE, ms=9, label="z")
        axes.set(xlim=(-3.2, 3.2), ylim=(-3.2, 3.2), aspect="equal", xlabel="Re z", ylabel="Im z")
        axes.set_title("z plane with the unit circle", fontsize=10, loc="left")
        axes.grid(True, color="#e1e0d9", lw=.8)
        axes.legend(fontsize=8, frameon=False, loc="upper right")
        self.cd_canvas.draw_idle()

    # -- rim explorer

    def build_rim(self, page):
        left, right = self.panel(page)
        ttk.Label(left, font=MONO, wraplength=380, justify="left",
                  text="N(K) = H·(-1 + (1-t)K) + (i/λ)(K - 5/2)²\nH = 1 - (1 + e^{iπ/4})·s").grid(
            row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
        change = lambda: self.later("rim", self.update_rim, 40)
        self.r = {"t": Control(left, 1, "t = y/x", .8, 1.2, 1., change),
                  "lam": Control(left, 2, "λ = x/(εA)", 20, 1400, 470, change, digits=2),
                  "s": Control(left, 3, "auxiliary s", 0, 1, .27, change)}
        buttons = ttk.Frame(left)
        buttons.grid(row=4, column=0, columnspan=3, sticky="w", pady=8)
        upper, lower = solver.collision_points()
        for column, (text, point) in enumerate((("Collision K = +50/3", upper),
                                                ("Collision K = -50/3", lower),
                                                ("One-meridian band", {"t": 1., "lam": 470., "s": .27}))):
            ttk.Button(buttons, text=text, command=lambda p=point: self.jump(p)).grid(
                row=0, column=column, padx=(0, 6))
        self.rim_table = Table(left, ("Quantity", "Value"), (150, 250), 8)
        self.rim_table.grid(row=5, column=0, columnspan=3, sticky="ew")
        self.rim_figure, self.rim_canvas = self.figure(right, 0, (5.2, 5.2))
        self.rim_note = self.note(right, 1)

    def jump(self, point):
        for key in ("t", "lam", "s"):
            self.r[key].set(point[key], notify=False)
        self.update_rim()

    def update_rim(self):
        try:
            r = solver.rim_slice(self.r["t"].get(), self.r["lam"].get(), self.r["s"].get())
        except ValueError as error:
            self.rim_note.configure(text=f"Solver error: {error}", foreground=BAD)
            return
        self.results["rim"] = r
        names = ("K+", "K-")
        rows = [(name, cplx(root), "") for name, root in zip(names, r["roots"])]
        rows += [(f"|{name}|·R", f"{modulus:.5f}", "") for name, modulus in zip(names, r["scaled_modulus"])]
        rows += [("Root separation", f"{r['separation']:.5f}", ""),
                 ("Roots inside rim", str(r["inside"]), "")]
        if r["single_periods"] is None:
            rows.append(("Single meridians", "divergent at the double root", "bad"))
        else:
            rows.append(("Single K+ / K-", " , ".join(f"{abs(complex(*v)):.3e}" for v in r["single_periods"]), ""))
        rows.append(("Paired period", cplx_e(r["paired_period"]), ""))
        self.rim_table.fill(rows)
        self.rim_note.configure(text=REGIMES[r["regime"]],
                                foreground=BAD if r["regime"] == "collision" else "")
        self.rim_figure.clear()
        axes = self.rim_figure.add_subplot()
        radius = r["radius"]
        axes.add_patch(_circle(1., GREY))
        axes.plot([0, radius, 2.5*radius], [0, 0, 0], "o", color=GREY, ms=4, label="punctures K = 0, 1, 5/2")
        for name, root, colour in zip(names, r["roots"], (BLUE, ORANGE)):
            x, y = root[0]*radius, root[1]*radius
            shrink = max(1., (x*x+y*y)**.5/1.6)
            axes.plot([x/shrink], [y/shrink], "o", color=colour, ms=10,
                      label=name+(" (off scale)" if shrink > 1 else ""))
        axes.set(xlim=(-1.75, 1.75), ylim=(-1.75, 1.75), aspect="equal",
                 xlabel="Re K · R", ylabel="Im K · R")
        axes.set_title("Base cap, rim |K| = 50/3 drawn as the unit circle", fontsize=10, loc="left")
        axes.grid(True, color="#e1e0d9", lw=.8)
        axes.legend(fontsize=8, frameon=False, loc="upper right")
        self.rim_canvas.draw_idle()


def _circle(radius, colour):
    from matplotlib.patches import Circle
    return Circle((0, 0), radius, fill=False, color=colour, lw=1.8)


def run() -> None:
    try:
        root = tk.Tk()
    except tk.TclError as error:
        raise SystemExit(f"No display is available for the desktop window ({error}). "
                         "Use the functional, edge, two-loop, contour, cd and rim "
                         "subcommands instead.") from None
    SolverApp(root)
    root.mainloop()


if __name__ == "__main__":
    run()
