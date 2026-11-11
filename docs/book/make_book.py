"""Build "Tiny book of EMT simulation" (docs/book/TinyBookOfEMT.pdf).

    python docs/book/make_book.py            figures + HTML + PDF (PDF printed with Microsoft Edge)
    python docs/book/make_book.py --no-pdf   figures + HTML only

Every number in the book is computed here by a small re-implementation of the EMTSim time
loop (the same companion models, the same order of operations), and checked against the
C++ solver (build/emtsim.exe) when it has been built.
"""
import math
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, Polygon, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
FIG = os.path.join(HERE, "fig")
NETS = os.path.join(HERE, "circuits")

INK = "#1a1a1a"
ACCENT = "#C55A11"      # simulation results and numbers (the orange-brown of the neural network book)
GREY = "#8a8a8a"
FAULT = "#f3dcdc"

plt.rcParams.update({
    "font.family": "Georgia", "font.size": 11, "svg.fonttype": "path",
    "axes.edgecolor": "#555", "axes.linewidth": 0.8, "axes.labelcolor": INK,
    "xtick.color": "#444", "ytick.color": "#444", "axes.spines.top": False, "axes.spines.right": False,
    "mathtext.fontset": "cm",
})


# ============================================================ a tiny EMT solver (same steps as EMTSim)
class Circuit:
    """Elements: ("R"|"L"|"C", a, b, value, name), ("V", node, peak, freq, name) ideal source to ground,
    ("SW", a, b, closed(t) -> bool, name) with Ron = 0.01, Roff = 1e6 as in EMTSim."""

    def __init__(self, n, elements, dt):
        self.n, self.el, self.dt = n, elements, dt
        self.known = [e[1] for e in elements if e[0] == "V"]
        self.unknown = [k for k in range(1, n + 1) if k not in self.known]
        self.state = {e[-1]: {"I": 0.0, "V": 0.0, "h": 0.0} for e in elements}

    def reff(self, e, t):
        k = e[0]
        if k == "R":
            return e[3]
        if k == "L":
            return 2 * e[3] / self.dt
        if k == "C":
            return self.dt / (2 * e[3])
        if k == "SW":
            return 0.01 if e[3](t) else 1.0e6

    def G(self, t):
        G = np.zeros((self.n + 1, self.n + 1))
        for e in self.el:
            if e[0] == "V":
                continue
            g = 1 / self.reff(e, t)
            a, b = e[1], e[2]
            G[a, a] += g; G[b, b] += g; G[a, b] -= g; G[b, a] -= g
        return G[1:, 1:]

    def step(self, t, euler=False):
        """One step at time t (a backward-Euler half step with euler=True: the history currents carry
        no voltage term, the conductances are the same). Returns dict with everything the book prints."""
        G = self.G(t)
        inj = np.zeros(self.n + 1)
        hist = {}
        for e in self.el:
            s = self.state[e[-1]]
            if e[0] == "L":
                s["h"] = s["I"] if euler else s["I"] + s["V"] / self.reff(e, t)
            elif e[0] == "C":
                s["h"] = -s["V"] / self.reff(e, t) if euler else -s["I"] - s["V"] / self.reff(e, t)
            else:
                s["h"] = 0.0
            if e[0] in "LC":
                hist[e[-1]] = s["h"]
                inj[e[1]] -= s["h"]; inj[e[2]] += s["h"]
        inj = inj[1:]
        U = [k - 1 for k in self.unknown]
        K = [k - 1 for k in self.known]
        vK = np.array([e[2] * math.cos(2 * math.pi * e[3] * t) for e in self.el if e[0] == "V"])
        rhs = inj[U] - G[np.ix_(U, K)] @ vK
        vU = np.linalg.solve(G[np.ix_(U, U)], rhs)
        v = np.zeros(self.n + 1)
        v[[k + 1 for k in K]] = vK
        v[[k + 1 for k in U]] = vU
        cur = {}
        for e in self.el:
            s = self.state[e[-1]]
            if e[0] == "V":
                continue
            s["V"] = v[e[1]] - v[e[2]]
            s["I"] = s["V"] / self.reff(e, t) + s["h"]
            cur[e[-1]] = s["I"]
        return {"t": t, "v": v, "i": cur, "hist": hist, "G": G, "rhs": rhs, "inj": inj, "vK": vK}

    def run(self, tstop):
        out, t = [], 0.0
        nsteps = int(round(tstop / self.dt))
        for k in range(1, nsteps + 1):
            out.append(self.step(k * self.dt))
        return out


def rl_circuit(dt, opens=None):
    if opens is None:
        return Circuit(2, [("V", 1, 10.0, 0.0, "BATTERY"), ("R", 1, 2, 1.0, "R"), ("L", 2, 0, 1e-3, "L")], dt)
    return Circuit(3, [("V", 1, 10.0, 0.0, "BATTERY"), ("SW", 1, 2, lambda t: t < opens - 1e-12, "SWITCH"),
                       ("R", 2, 3, 1.0, "R"), ("L", 3, 0, 1e-3, "L")], dt)


def rl_open_cda(dt, t_open=3e-3, tstop=5e-3):
    """The circuit of chapter 8 with chatter removal: the switch opens at t_open (a time point) and the
    next two steps are backward-Euler half steps. Returns the grid steps and the half-step solution."""
    c = Circuit(3, [("V", 1, 10.0, 0.0, "BATTERY"), ("SW", 1, 2, lambda t: t <= t_open + 1e-12, "SWITCH"),
                    ("R", 2, 3, 1.0, "R"), ("L", 3, 0, 1e-3, "L")], dt)
    out, half = [], None
    k, n = 1, int(round(tstop / dt))
    while k <= n:
        t = k * dt
        out.append(c.step(t))
        k += 1
        if abs(t - t_open) < 1e-12:
            half = c.step(t + 0.5 * dt, euler=True)
            out.append(c.step(t + dt, euler=True))
            k += 1
    return out, half


def network(dt, fault=None):
    el = [("V", 1, 100.0, 50.0, "SOURCE"), ("R", 1, 2, 1.0, "R1"), ("L", 2, 3, 0.01, "L"),
          ("C", 3, 0, 10e-6, "C"), ("R", 3, 0, 100.0, "LOAD")]
    if fault:
        el.append(("SW", 3, 0, lambda t: fault[0] - 1e-12 <= t < fault[1] - 1e-12, "FAULT"))
    return Circuit(3, el, dt)


# ------------------------------------------------------------ check against the C++ solver
def emtsim(netfile, *flags):
    exe = os.path.join(ROOT, "build", "emtsim.exe" if os.name == "nt" else "emtsim")
    if not os.path.exists(exe):
        return None
    d = tempfile.mkdtemp()
    try:
        subprocess.run([exe, netfile, *flags], cwd=d, check=True, capture_output=True)
        v = np.loadtxt(os.path.join(d, "results.csv"), delimiter=",", skiprows=1)
        with open(os.path.join(d, "currents.csv")) as f:
            names = f.readline().strip().split(",")
        c = np.loadtxt(os.path.join(d, "currents.csv"), delimiter=",", skiprows=1)
        return v, dict(zip(names, c.T))
    finally:
        shutil.rmtree(d, ignore_errors=True)


def check(name, steps, current_names, *flags):
    res = emtsim(os.path.join(NETS, name), *flags)
    if res is None:
        print("  (build/emtsim.exe not found: %s not checked)" % name)
        return
    v, c = res
    mine_v = np.array([s["v"][1:] for s in steps])
    dv = np.abs(mine_v - v[:, 1:]).max()
    di = max(np.abs(np.array([s["i"][n] for s in steps]) - c[n]).max() for n in current_names)
    scale = max(1.0, np.abs(v[:, 1:]).max())
    ok = dv < 1e-9 * scale and di < 1e-9 * scale
    print("  %-18s EMTSim agrees: max difference %.1e V, %.1e A  %s" % (name, dv, di, "OK" if ok else "MISMATCH"))
    if not ok:
        raise SystemExit("book numbers do not match EMTSim for " + name)


# ============================================================ schematic drawing
LW = 1.6


def _frame(p1, p2, body):
    p1, p2 = np.asarray(p1, float), np.asarray(p2, float)
    d = p2 - p1
    L = np.hypot(*d)
    u = d / L
    n = np.array([-u[1], u[0]])
    a = p1 + u * (L - body) / 2
    T = lambda s, w: a + u * s + n * w
    return p1, p2, a, a + u * body, u, n, T


def line(ax, *pts, **kw):
    pts = np.asarray(pts, float)
    ax.plot(pts[:, 0], pts[:, 1], color=kw.get("color", INK), lw=kw.get("lw", LW), solid_capstyle="round",
            zorder=kw.get("zorder", 2))


def resistor(ax, p1, p2, body=1.2):
    p1, p2, a, b, u, n, T = _frame(p1, p2, body)
    s = np.linspace(0, body, 14)
    w = np.zeros_like(s)
    w[1:-1] = [0.17 if k % 2 else -0.17 for k in range(1, len(s) - 1)]
    line(ax, p1, a); line(ax, b, p2)
    line(ax, *[T(si, wi) for si, wi in zip(s, w)])


def inductor(ax, p1, p2, body=1.2, side=1):
    p1, p2, a, b, u, n, T = _frame(p1, p2, body)
    line(ax, p1, a); line(ax, b, p2)
    r = body / 8
    pts = []
    for k in range(4):
        c = (2 * k + 1) * r
        for th in np.linspace(math.pi, 0, 16):
            pts.append(T(c + r * math.cos(th), side * r * 1.3 * math.sin(th)))
    line(ax, *pts)


def capacitor(ax, p1, p2, gap=0.14, half=0.32):
    p1, p2, a, b, u, n, T = _frame(p1, p2, 0)
    m = (p1 + p2) / 2
    line(ax, p1, m - u * gap); line(ax, m + u * gap, p2)
    line(ax, m - u * gap - n * half, m - u * gap + n * half, lw=2.2)
    line(ax, m + u * gap - n * half, m + u * gap + n * half, lw=2.2)


def vsource(ax, p1, p2, kind="dc", r=0.42):
    """p1 = minus terminal, p2 = plus terminal"""
    p1, p2, a, b, u, n, T = _frame(p1, p2, 2 * r)
    c = (p1 + p2) / 2
    line(ax, p1, c - u * r); line(ax, c + u * r, p2)
    ax.add_patch(Circle(c, r, fill=False, ec=INK, lw=LW, zorder=2))
    if kind == "ac":
        s = np.linspace(-0.6, 0.6, 40) * r
        pts = [c + n * si + u * 0.35 * r * math.sin(si / r / 0.6 * math.pi) for si in s]
        line(ax, *pts, lw=1.3)
    else:
        ax.text(*(c + u * r * 0.45), "+", ha="center", va="center", fontsize=11, color=INK)
        ax.text(*(c - u * r * 0.5), "−", ha="center", va="center", fontsize=11, color=INK)


def isource(ax, p1, p2, r=0.36):
    """arrow from p1 to p2"""
    p1, p2, a, b, u, n, T = _frame(p1, p2, 2 * r)
    c = (p1 + p2) / 2
    line(ax, p1, c - u * r); line(ax, c + u * r, p2)
    ax.add_patch(Circle(c, r, fill=False, ec=INK, lw=LW, zorder=2))
    ax.annotate("", xy=c + u * r * 0.7, xytext=c - u * r * 0.7,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.4, mutation_scale=11))


def switch(ax, p1, p2, closed=False, body=1.0):
    p1, p2, a, b, u, n, T = _frame(p1, p2, body)
    line(ax, p1, a); line(ax, b, p2)
    for q in (a, b):
        ax.add_patch(Circle(q, 0.07, fc="white", ec=INK, lw=1.4, zorder=3))
    tip = b if closed else a + u * body * 0.92 + n * 0.45
    line(ax, a, tip)


def ground(ax, p, down=0.35):
    p = np.asarray(p, float)
    line(ax, p, p - [0, down])
    for k, w in enumerate((0.32, 0.2, 0.08)):
        y = p[1] - down - k * 0.1
        line(ax, (p[0] - w, y), (p[0] + w, y), lw=1.4)


def node(ax, p, label=None, dx=0.0, dy=0.32):
    ax.add_patch(Circle(p, 0.07, fc=INK, ec=INK, zorder=4))
    if label is not None:
        q = (p[0] + dx, p[1] + dy)
        ax.add_patch(Circle(q, 0.17, fc=INK, ec=INK, zorder=4))
        ax.text(*q, str(label), color="white", ha="center", va="center", fontsize=8.5, zorder=5, fontweight="bold")


def text(ax, p, s, **kw):
    kw.setdefault("ha", "center"); kw.setdefault("va", "center")
    kw.setdefault("fontsize", 11); kw.setdefault("color", INK)
    ax.text(p[0], p[1], s, **kw)


def canvas(w, h, xlim, ylim):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal"); ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1)
    return fig, ax


def save(fig, name):
    fig.savefig(os.path.join(FIG, name), transparent=False, facecolor="white")
    plt.close(fig)


# ============================================================ figures
def fig_cover():
    """Line art for the cover tile: a steady sine wave, a switching transient, and the time steps."""
    fig, ax = plt.subplots(figsize=(2.3, 5.4))
    fig.subplots_adjust(0, 0, 1, 1)
    ax.axis("off")
    t = np.linspace(0, 1, 2000)
    # vertical "time" axis going down the tile, waveform swinging left-right
    y = 1 - t
    env = np.where(t < 0.38, 0, np.exp(-(t - 0.38) * 7))
    x = 0.55 * np.sin(2 * np.pi * 3 * t) * (1 - 0.6 * (t >= 0.38)) + 0.9 * env * np.sin(2 * np.pi * 26 * (t - 0.38))
    ax.plot(x, y, color=INK, lw=1.5)
    ts = t[::60]
    xs = np.interp(ts, t, x)
    ax.scatter(xs, 1 - ts, s=7, color=INK, zorder=3)
    for k in np.linspace(0, 1, 34):
        ax.plot([-1.25, -1.1], [k, k], color=INK, lw=0.6)
    ax.plot([-1.175, -1.175], [0, 1], color=INK, lw=0.6)
    ax.set_xlim(-1.35, 1.25); ax.set_ylim(-0.04, 1.04)
    save(fig, "cover.svg")


def fig_photo_movie():
    fig = plt.figure(figsize=(6.6, 2.5))
    # the photo: one phasor
    ax = fig.add_axes([0.0, 0.08, 0.26, 0.84]); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(FancyBboxPatch((-1.2, -1.2), 2.4, 2.4, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fill=False, ec=INK, lw=1.4))
    ax.add_patch(Circle((0, 0), 0.85, fill=False, ec=GREY, lw=0.8, ls="--"))
    ax.annotate("", xy=(0.6, 0.6), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=ACCENT, lw=2, mutation_scale=14))
    ax.text(0, -1.5, "Phasor: one photo", ha="center", fontsize=10.5, fontweight="bold")
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.75, 1.3)
    # the movie: a film strip
    ax = fig.add_axes([0.3, 0.08, 0.7, 0.84]); ax.axis("off")
    nf = 6
    W, H = 1.0, 0.8
    ax.add_patch(plt.Rectangle((-0.1, -0.25), nf * (W + 0.1) + 0.1, H + 0.5, fc=INK, ec=INK))
    for k in range(nf * 3 + 1):
        ax.add_patch(plt.Rectangle((-0.06 + k * (W + 0.1) / 3, H + 0.1), 0.12, 0.08, fc="white", ec="white"))
        ax.add_patch(plt.Rectangle((-0.06 + k * (W + 0.1) / 3, -0.18), 0.12, 0.08, fc="white", ec="white"))
    tt = np.linspace(0, 1, 300)
    wave = lambda x: np.where(x < 0.3, 0.0, np.sin(2 * np.pi * 2.2 * (x - 0.3)) * (1 - np.exp(-(x - 0.3) * 9))
                              + 0.6 * np.exp(-(x - 0.3) * 7) * np.sin(2 * np.pi * 9 * (x - 0.3)))
    for k in range(nf):
        x0 = k * (W + 0.1)
        ax.add_patch(plt.Rectangle((x0, 0), W, H, fc="white", ec="white"))
        tk = (k + 1) / nf
        m = tt <= tk
        ax.plot(x0 + 0.05 + 0.9 * tt[m], H / 2 + 0.3 * wave(tt[m]), color=GREY, lw=0.9)
        ax.plot(x0 + 0.05 + 0.9 * tk, H / 2 + 0.3 * wave(np.array([tk]))[0], "o", color=ACCENT, ms=4.5)
        ax.text(x0 + W / 2, -0.47, ["t = Δt", "2Δt", "3Δt", "4Δt", "5Δt", "6Δt"][k], ha="center", fontsize=9.5)
    ax.text(nf * (W + 0.1) / 2, -0.82, "EMT simulation: a movie, one frame every Δt", ha="center",
            fontsize=10.5, fontweight="bold")
    ax.set_xlim(-0.15, nf * (W + 0.1) + 0.05); ax.set_ylim(-1.0, H + 0.3)
    save(fig, "photo_movie.svg")


def fig_rl():
    fig, ax = canvas(5.6, 2.9, (-1.9, 9.4), (-1.2, 4.4))
    vsource(ax, (0, 0), (0, 3))
    switch(ax, (0, 3), (3, 3), closed=False)
    resistor(ax, (3, 3), (6, 3))
    line(ax, (6, 3), (8, 3))
    inductor(ax, (8, 3), (8, 0))
    line(ax, (0, 0), (8, 0))
    ground(ax, (4, 0))
    node(ax, (0, 3), 1, dy=0.35); node(ax, (6, 3), 2, dy=0.35)
    text(ax, (-0.55, 1.5), "V\n10 V", ha="right", fontsize=10)
    text(ax, (1.5, 3.75), "switch closes\nat t = 0", fontsize=9.5)
    text(ax, (4.5, 3.6), "R = 1 Ω", fontsize=10)
    text(ax, (8.5, 1.5), "L\n1 mH", ha="left", fontsize=10)
    ax.annotate("", xy=(7.6, 3.25), xytext=(6.8, 3.25), arrowprops=dict(arrowstyle="-|>", color=ACCENT, lw=1.4))
    text(ax, (7.2, 3.6), "i", color=ACCENT, fontsize=12, style="italic")
    save(fig, "rl_circuit.svg")


def fig_trapezoid():
    fig, ax = plt.subplots(figsize=(3.4, 2.3))
    t = np.linspace(0, 1, 200)
    v = 0.3 + 0.75 * t + 0.12 * np.sin(5 * t)
    ax.plot(t, v, color=INK, lw=1.6)
    a, b = 0.25, 0.75
    va, vb = np.interp([a, b], t, v)
    ax.add_patch(Polygon([(a, 0), (a, va), (b, vb), (b, 0)], fc=ACCENT, alpha=0.22, ec=ACCENT, lw=1.2))
    m = (t >= a) & (t <= b)
    ax.fill_between(t[m], np.interp(t[m], [a, b], [va, vb]), v[m], color=GREY, alpha=0.3, lw=0)
    ax.plot([a, b], [va, vb], color=ACCENT, lw=1.4)
    ax.set_xticks([a, b]); ax.set_xticklabels(["t − Δt", "t"])
    ax.set_yticks([va, vb]); ax.set_yticklabels(["v(t − Δt)", "v(t)"])
    ax.set_ylim(0, 1.25); ax.set_xlim(0, 1)
    ax.text((a + b) / 2, 0.2, "area ≈ Δt × (average\nof the two heights)", ha="center", fontsize=9)
    fig.tight_layout()
    save(fig, "trapezoid.svg")


def fig_companion():
    fig, ax = canvas(6.4, 4.1, (-0.6, 13.4), (-0.9, 7.6))
    # inductor
    for y0, kind in ((4.3, "L"), (0.0, "C")):
        node(ax, (0, y0 + 2.6)); node(ax, (0, y0))
        text(ax, (0, y0 + 2.95), "k", style="italic"); text(ax, (0, y0 - 0.38), "m", style="italic")
        if kind == "L":
            inductor(ax, (0, y0 + 2.6), (0, y0))
            text(ax, (0.55, y0 + 1.3), "L", ha="left", style="italic", fontsize=12)
        else:
            capacitor(ax, (0, y0 + 2.6), (0, y0))
            text(ax, (0.55, y0 + 1.3), "C", ha="left", style="italic", fontsize=12)
        ax.annotate("", xy=(3.6, y0 + 1.3), xytext=(1.9, y0 + 1.3),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6, mutation_scale=16))
        text(ax, (2.75, y0 + 1.75), "one time step", fontsize=8.5, color="#444")
        x1, x2 = 5.2, 7.6
        node(ax, (6.4, y0 + 2.6)); node(ax, (6.4, y0))
        text(ax, (6.4, y0 + 2.95), "k", style="italic"); text(ax, (6.4, y0 - 0.38), "m", style="italic")
        line(ax, (x1, y0 + 2.6), (x2, y0 + 2.6)); line(ax, (x1, y0), (x2, y0))
        resistor(ax, (x1, y0 + 2.6), (x1, y0), body=1.1)
        isource(ax, (x2, y0 + 2.6), (x2, y0))
        if kind == "L":
            text(ax, (x2 + 0.6, y0 + 1.75), r"$R_L = 2L\,/\,\Delta t$", ha="left", fontsize=11)
            text(ax, (x2 + 0.6, y0 + 0.85), r"$I_{\rm hist}(t) = i(t-\Delta t) + v(t-\Delta t)\,/\,R_L$", ha="left", fontsize=11)
        else:
            text(ax, (x2 + 0.6, y0 + 1.75), r"$R_C = \Delta t\,/\,2C$", ha="left", fontsize=11)
            text(ax, (x2 + 0.6, y0 + 0.85), r"$I_{\rm hist}(t) = -i(t-\Delta t) - v(t-\Delta t)\,/\,R_C$", ha="left", fontsize=11)
        text(ax, (x1 + 0.45, y0 + 1.3), r"$R_%s$" % kind, ha="left", fontsize=11)
    save(fig, "companion.svg")


def fig_rl_companion():
    fig, ax = canvas(5.8, 2.6, (-1.9, 10.2), (-1.2, 3.8))
    vsource(ax, (0, 0), (0, 3))
    resistor(ax, (0, 3), (4, 3))
    line(ax, (4, 3), (8.4, 3))
    resistor(ax, (6.2, 3), (6.2, 0), body=1.1)
    isource(ax, (8.4, 3), (8.4, 0))
    line(ax, (0, 0), (8.4, 0))
    ground(ax, (3, 0))
    node(ax, (0, 3), 1, dy=0.35); node(ax, (4, 3), 2, dy=0.35)
    text(ax, (-0.55, 1.5), "V\n10 V", ha="right", fontsize=10)
    text(ax, (2, 3.55), "R = 1 Ω", fontsize=10)
    text(ax, (6.65, 1.5), "$R_L$\n20 Ω", ha="left", fontsize=10)
    text(ax, (8.95, 1.5), r"$I_{\rm hist}$", ha="left", fontsize=11)
    ax.add_patch(FancyBboxPatch((5.35, -0.45), 4.55, 3.9, boxstyle="round,pad=0,rounding_size=0.25",
                                fill=False, ec=ACCENT, lw=1.1, ls="--"))
    text(ax, (7.6, -0.8), "the inductor, for one time step", fontsize=9, color=ACCENT)
    save(fig, "rl_companion.svg")


def fig_network():
    fig, ax = canvas(6.4, 2.7, (-2.4, 13.2), (-1.3, 4.3))
    vsource(ax, (0, 0), (0, 3), kind="ac")
    resistor(ax, (0, 3), (3.4, 3))
    inductor(ax, (3.4, 3), (7.2, 3))
    line(ax, (7.2, 3), (11.5, 3))
    capacitor(ax, (8.8, 3), (8.8, 0))
    resistor(ax, (11.5, 3), (11.5, 0))
    line(ax, (0, 0), (11.5, 0))
    ground(ax, (5, 0))
    node(ax, (0, 3), 1, dy=0.38); node(ax, (3.4, 3), 2, dy=0.38); node(ax, (7.2, 3), 3, dy=0.38)
    node(ax, (8.8, 3)); node(ax, (8.8, 0)); node(ax, (7.2, 0))
    text(ax, (-0.6, 1.5), "100 V peak\n50 Hz", ha="right", fontsize=9.5)
    text(ax, (1.7, 3.6), "R1 = 1 Ω", fontsize=10)
    text(ax, (5.3, 3.75), "L = 10 mH", fontsize=10)
    text(ax, (9.25, 1.5), "C\n10 μF", ha="left", fontsize=10)
    text(ax, (11.95, 1.5), "LOAD\n100 Ω", ha="left", fontsize=10)
    # fault (drawn dashed in chapter 5 only)
    save(fig, "network.svg")
    # the same with the fault
    fig, ax = canvas(6.4, 2.7, (-2.4, 13.2), (-1.3, 4.3))
    vsource(ax, (0, 0), (0, 3), kind="ac")
    resistor(ax, (0, 3), (3.4, 3))
    inductor(ax, (3.4, 3), (7.2, 3))
    line(ax, (7.2, 3), (11.5, 3))
    capacitor(ax, (8.8, 3), (8.8, 0))
    resistor(ax, (11.5, 3), (11.5, 0))
    line(ax, (0, 0), (11.5, 0))
    ground(ax, (5, 0))
    switch(ax, (7.2, 3), (7.2, 0), body=1.0)
    node(ax, (0, 3), 1, dy=0.38); node(ax, (3.4, 3), 2, dy=0.38); node(ax, (7.2, 3), 3, dy=0.38)
    text(ax, (6.95, 1.5), "FAULT", ha="right", fontsize=9.5, color="#B03030", fontweight="bold")
    text(ax, (-0.6, 1.5), "100 V peak\n50 Hz", ha="right", fontsize=9.5)
    text(ax, (1.7, 3.6), "R1 = 1 Ω", fontsize=10)
    text(ax, (5.3, 3.75), "L = 10 mH", fontsize=10)
    text(ax, (9.25, 1.5), "C\n10 μF", ha="left", fontsize=10)
    text(ax, (11.95, 1.5), "LOAD\n100 Ω", ha="left", fontsize=10)
    save(fig, "network_fault.svg")


def fig_loop():
    """What is known and what is unknown in one time step."""
    fig, ax = canvas(7.4, 3.2, (-0.2, 16.0), (-1.9, 5.6))
    def box(x, y, w, h, title, body, fc="white", ec=INK):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.25", fc=fc, ec=ec, lw=1.4))
        ax.text(x + w / 2, y + h - 0.45, title, ha="center", va="center", fontsize=10.5, fontweight="bold")
        ax.text(x + w / 2, y + h / 2 - 0.35, body, ha="center", va="center", fontsize=9.5, linespacing=1.35)
    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.4, mutation_scale=14))
    box(0.0, 2.2, 3.6, 2.6, "Last step (known)", "v(t − Δt)\ni(t − Δt)", fc="#F3F4F6")
    box(4.8, 2.2, 5.4, 2.6, "The notes", "I_hist(t) = i(t−Δt) + v(t−Δt) / R_L\nfor every L and C", fc=NOTE, ec="#C9B44A")
    box(11.4, 2.2, 4.4, 2.6, "Solve", "G v(t) = i_hist(t)\nunknown: v(t)")
    box(11.4, -0.4, 4.4, 1.9, "Branch currents", "i(t) = v(t) / R_eq + I_hist(t)")
    arrow(3.6, 3.5, 4.8, 3.5)
    arrow(10.2, 3.5, 11.4, 3.5)
    arrow(13.6, 2.2, 13.6, 1.5)
    line(ax, (13.6, -0.4), (13.6, -1.1), (1.8, -1.1), lw=1.4)
    arrow(1.8, -1.1, 1.8, 2.2)
    ax.text(7.9, -1.45, "this step's voltages and currents become the next step's last step", ha="center", va="center",
            fontsize=9, color=ACCENT)
    ax.text(7.9, 5.3, "G never changes; only the right-hand side does", ha="center", fontsize=9.5, style="italic", color="#444")
    save(fig, "loop.svg")


def fig_rl_open():
    fig, ax = canvas(5.6, 2.9, (-1.9, 9.4), (-1.2, 4.4))
    vsource(ax, (0, 0), (0, 3))
    switch(ax, (0, 3), (3, 3), closed=True)
    resistor(ax, (3, 3), (6, 3))
    line(ax, (6, 3), (8, 3))
    inductor(ax, (8, 3), (8, 0))
    line(ax, (0, 0), (8, 0))
    ground(ax, (4, 0))
    text(ax, (-0.55, 1.5), "V\n10 V", ha="right", fontsize=10)
    text(ax, (1.5, 3.75), "switch opens\nat t = 3 ms", fontsize=9.5)
    text(ax, (4.5, 3.6), "R = 1 Ω", fontsize=10)
    text(ax, (8.5, 1.5), "L\n1 mH", ha="left", fontsize=10)
    save(fig, "rl_open_circuit.svg")


# ============================================================ HTML helpers (MathML)
import re

_TOK = re.compile(r"&#\d+;|'[A-Za-zΩμ]+|\d+(?:\.\d+)?|[A-Za-zΔτσπω]+|\S")


def X(s):
    """plain text -> MathML tokens; strings that already are MathML pass through"""
    if s.startswith("<"):
        return s
    out = []
    for tok in _TOK.findall(s):
        if tok.startswith("&#"):
            out.append("<mspace width='1em'/>")
        elif tok.startswith("'"):          # unit: upright, thin space before it
            out.append("<mspace width='0.2em'/><mi mathvariant='normal'>%s</mi>" % tok[1:])
        elif tok in "()":
            out.append("<mo stretchy='false'>%s</mo>" % tok)
        elif re.fullmatch(r"\d+(?:\.\d+)?", tok):
            out.append("<mn>%s</mn>" % tok)
        elif re.fullmatch(r"[A-Za-zΔτσπω]+", tok):
            if tok in ("Δt",):
                out.append("<mi mathvariant='normal'>Δ</mi><mi>t</mi>")
            elif len(tok) == 1 or tok == "Δ":
                out.append("<mi>%s</mi>" % tok)
            else:
                out.append("<mi mathvariant='normal'>%s</mi>" % tok)
        else:
            out.append("<mo>%s</mo>" % tok)
    return "".join(out)


def cat(*parts):
    return "<mrow>" + "".join(X(p) for p in parts) + "</mrow>"


def fr(a, b):
    return "<mfrac>%s%s</mfrac>" % (cat(a), cat(b))


def sb(a, b):
    return "<msub>%s%s</msub>" % (cat(a), cat(b))


def sp(a, b):
    return "<msup>%s%s</msup>" % (cat(a), cat(b))


def br(*parts, o="[", c="]"):
    return "<mrow><mo>%s</mo>%s<mo>%s</mo></mrow>" % (o, cat(*parts), c)


def integral(lo, hi, *body):
    return "<mrow><msubsup><mo>∫</mo>%s%s</msubsup>%s</mrow>" % (cat(lo), cat(hi), cat(*body))


def num(x, d=4):
    if isinstance(x, str):
        return "<mn class='num'>%s</mn>" % x
    s = ("%.*f" % (d, x)).replace("-", "−")
    if s in ("−0." + "0" * d,):
        s = s[1:]
    return "<mn class='num'>%s</mn>" % s


def mat(A, d=4):
    A = np.atleast_2d(A)
    rows = "".join("<mtr>" + "".join("<mtd>%s</mtd>" % num(x, d) for x in r) + "</mtr>" for r in A)
    return "<mrow><mo>[</mo><mtable columnspacing='1.1em'>%s</mtable><mo>]</mo></mrow>" % rows


def col(v, d=4):
    return mat(np.asarray(v).reshape(-1, 1), d)


def eq(*parts, label=None):
    s = "<math display='block'>%s</math>" % cat(*parts)
    if label:
        return "<div class='eq'>%s<span class='eqn'>(%s)</span></div>" % (s, label)
    return "<div class='eq'>%s</div>" % s


def m(*parts):
    """inline maths"""
    return "<math>%s</math>" % cat(*parts)


I_HIST = sb("I", "hist")
BG, BV, BI = ("<mi mathvariant='normal'>%s</mi>" % x for x in ("𝐆", "𝐯", "𝐢"))   # bold G, v, i
R_L = sb("R", "L")
R_C = sb("R", "C")


def f(x, d=4):
    return ("%.*f" % (d, x)).replace("-", "−")


# ============================================================ cartoons
from matplotlib.patches import Ellipse, Rectangle
from matplotlib.path import Path
from matplotlib.patches import PathPatch

HAND = "Segoe Print"
NOTE = "#FFF3A6"


def stick(ax, x, y, s=1.0, arms="down", look=1):
    """a stick figure like the ones in the neural network book; look = +1 right, -1 left"""
    lw = 1.6
    hy = y + 2.15 * s
    ax.add_patch(Circle((x, hy), 0.36 * s, fc="white", ec=INK, lw=lw, zorder=5))
    for dx in (-0.12, 0.12):
        ax.add_patch(Circle((x + (dx + 0.05 * look) * s, hy + 0.06 * s), 0.035 * s, fc=INK, ec=INK, zorder=6))
    line(ax, (x, hy - 0.36 * s), (x, y + 0.9 * s), lw=lw, zorder=4)
    line(ax, (x, y + 0.9 * s), (x - 0.32 * s, y), lw=lw, zorder=4)
    line(ax, (x, y + 0.9 * s), (x + 0.32 * s, y), lw=lw, zorder=4)
    sh = (x, y + 1.55 * s)
    if arms == "down":
        line(ax, sh, (x - 0.45 * s, y + 1.0 * s), lw=lw); line(ax, sh, (x + 0.45 * s, y + 1.0 * s), lw=lw)
    elif arms == "hold":      # both hands forward, in the look direction
        line(ax, sh, (x + 0.55 * s * look, y + 1.3 * s), lw=lw)
        line(ax, sh, (x + 0.5 * s * look, y + 1.5 * s), lw=lw)
        line(ax, sh, (x - 0.4 * s * look, y + 1.0 * s), lw=lw)
    elif arms == "wave":
        line(ax, sh, (x + 0.5 * s * look, y + 2.35 * s), lw=lw)
        line(ax, sh, (x - 0.45 * s * look, y + 1.0 * s), lw=lw)
    elif arms == "think":
        line(ax, sh, (x + 0.3 * s * look, y + 1.25 * s), lw=lw)
        line(ax, (x + 0.3 * s * look, y + 1.25 * s), (x + 0.18 * s * look, y + 1.8 * s), lw=lw)
        line(ax, sh, (x - 0.45 * s * look, y + 1.0 * s), lw=lw)
    return (x + 0.55 * s * look, y + 1.35 * s)    # where a held object goes


def bubble(ax, cx, cy, w, h, txt, tail, size=9.5):
    """speech bubble centred at (cx, cy) with its tail pointing at `tail`"""
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h, boxstyle="round,pad=0.05,rounding_size=0.35",
                                fc="white", ec=INK, lw=1.2, zorder=7))
    tx, ty = tail
    bx = min(max(tx, cx - w / 2 + 0.4), cx + w / 2 - 0.4)
    by = cy - h / 2 if ty < cy else cy + h / 2
    ax.add_patch(Polygon([(bx - 0.22, by), (bx + 0.22, by), (tx, ty)], closed=True, fc="white", ec=INK, lw=1.2, zorder=6))
    ax.plot([bx - 0.2, bx + 0.2], [by, by], color="white", lw=2.6, zorder=8)
    ax.text(cx, cy, txt, ha="center", va="center", fontsize=size, family=HAND, color=INK, zorder=9, linespacing=1.15)


def eyes(ax, x, y, s=1.0, look=0.0):
    for dx in (-0.17, 0.17):
        ax.add_patch(Ellipse((x + dx * s, y), 0.24 * s, 0.3 * s, fc="white", ec=INK, lw=1.2, zorder=6))
        ax.add_patch(Circle((x + (dx + 0.04 * look) * s, y - 0.02 * s), 0.05 * s, fc=INK, zorder=7))


def smile(ax, x, y, s=1.0, sad=False):
    th = np.linspace(math.pi * 1.15, math.pi * 1.85, 20) if not sad else np.linspace(math.pi * 0.15, math.pi * 0.85, 20)
    r = 0.22 * s
    yy = y + (0 if not sad else -0.3 * s)
    ax.plot(x + r * np.cos(th), yy + r * np.sin(th) + (0.0 if not sad else 0.0), color=INK, lw=1.3, zorder=7)


def fig_cast():
    fig, ax = canvas(6.6, 3.0, (-0.2, 13.4), (-1.5, 4.6))
    cols = [(2.0, "R", "I live in the moment.\nNo memory at all!", "Resistor", r"$v = R\,i$"),
            (6.6, "L", "Change my current?\nNot so fast!", "Inductor", r"$v = L\,\dfrac{di}{dt}$"),
            (11.2, "C", "Change my voltage?\nNot so fast!", "Capacitor", r"$i = C\,\dfrac{dv}{dt}$")]
    for cx, kind, says, name, law in cols:
        if kind == "R":
            resistor(ax, (cx - 1.5, 0.6), (cx + 1.5, 0.6), body=1.6)
        elif kind == "L":
            inductor(ax, (cx - 1.5, 0.6), (cx + 1.5, 0.6), body=1.8)
        else:
            capacitor(ax, (cx - 1.5, 0.6), (cx + 1.5, 0.6), gap=0.16, half=0.42)
        eyes(ax, cx, 1.35, look=0.5)
        smile(ax, cx, 1.22, 0.8)
        bubble(ax, cx, 3.35, 3.9, 1.35, says, (cx + 0.2, 1.75), size=9.2)
        ax.text(cx, -0.25, name, ha="center", fontsize=11, fontweight="bold", color=INK)
        ax.text(cx, -1.0, law, ha="center", fontsize=12, color=ACCENT)
    save(fig, "cast.svg")


def fig_pete_emma():
    fig, ax = canvas(6.6, 3.0, (-0.5, 13.5), (-1.0, 5.2))
    # Phasor Pete with a photo
    h = stick(ax, 2.0, 0, 1.25, arms="hold", look=1)
    ax.add_patch(Rectangle((h[0] - 0.05, h[1] - 0.45), 1.3, 1.0, fc="white", ec=INK, lw=1.4, zorder=6))
    xs = np.linspace(0, 1.1, 60)
    ax.plot(h[0] + 0.05 + xs, h[1] + 0.05 + 0.3 * np.sin(2 * np.pi * 1.8 * xs), color=ACCENT, lw=1.3, zorder=7)
    bubble(ax, 2.6, 4.3, 4.4, 1.35, "A steady sine wave:\n100 V, 50 Hz. Done!", (2.3, 3.2))
    ax.text(2.0, -0.75, "Phasor Pete", ha="center", fontsize=10.5, fontweight="bold")
    # EMT Emma with a film strip
    h = stick(ax, 11.0, 0, 1.25, arms="hold", look=-1)
    x0 = h[0] - 2.3
    ax.add_patch(Rectangle((x0, h[1] - 0.45), 2.3, 0.95, fc=INK, ec=INK, zorder=6))
    for k in range(4):
        ax.add_patch(Rectangle((x0 + 0.1 + k * 0.55, h[1] - 0.3), 0.45, 0.65, fc="white", ec="white", zorder=7))
        tt = np.linspace(0, 1, 20)
        ax.plot(x0 + 0.12 + k * 0.55 + 0.4 * tt, h[1] + 0.03 + 0.22 * np.sin(2 * np.pi * (tt + k * 0.6)) * (0.4 + 0.5 * k / 3),
                color=ACCENT, lw=0.9, zorder=8)
    bubble(ax, 10.2, 4.3, 5.0, 1.35, "And what happens right\nafter the switch closes?", (10.8, 3.2))
    ax.text(11.0, -0.75, "EMT Emma", ha="center", fontsize=10.5, fontweight="bold")
    save(fig, "pete_emma.svg")


def fig_sticky(i1, v1):
    fig, ax = canvas(5.6, 2.5, (-0.5, 11.5), (-0.6, 4.6))
    inductor(ax, (0.3, 1.0), (4.3, 1.0), body=2.2)
    eyes(ax, 2.3, 1.85, look=1)
    smile(ax, 2.3, 1.7, 0.8)
    # the sticky note
    ax.add_patch(Polygon([(5.2, 0.0), (8.6, -0.15), (8.75, 3.0), (5.3, 3.15)], closed=True, fc=NOTE, ec="#C9B44A", lw=1, zorder=4))
    ax.text(6.98, 1.5, "Note to self\n(last step):\ni = %.4f A\nv = %.4f V" % (i1, v1), ha="center", va="center",
            family=HAND, fontsize=11.5, rotation=2, zorder=5)
    ax.annotate("", xy=(5.1, 1.2), xytext=(3.4, 1.9), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2,
                                                                      connectionstyle="arc3,rad=-0.3"))
    bubble(ax, 2.2, 3.75, 4.3, 1.2, "My memory: everything\nI need from the past.", (2.5, 2.3), size=10.5)
    save(fig, "sticky.svg")


def fig_stamp(gl):
    fig, ax = canvas(6.4, 2.7, (-0.3, 14.2), (-0.6, 5.6))
    # the stamp
    ax.add_patch(FancyBboxPatch((1.1, 4.0), 1.1, 1.2, boxstyle="round,pad=0.05,rounding_size=0.3", fc="#7a4a2a", ec=INK, lw=1.2))
    ax.add_patch(Rectangle((0.2, 3.5), 2.9, 0.5, fc="#a0673f", ec=INK, lw=1.2))
    ax.add_patch(Rectangle((0.2, 1.1), 2.9, 2.4, fc="white", ec=INK, lw=1.4))
    for (r, c), t in {(0, 0): "+g", (0, 1): "−g", (1, 0): "−g", (1, 1): "+g"}.items():
        ax.text(0.9 + c * 1.45, 2.85 - r * 1.15, t, ha="center", va="center", fontsize=13, color=ACCENT, fontweight="bold")
    ax.text(0.9, 0.7, "k", ha="center", style="italic"); ax.text(2.35, 0.7, "m", ha="center", style="italic")
    ax.text(1.65, 0.1, "the stamp of any\nconductance g", ha="center", va="top", fontsize=8.5, color="#444")
    ax.annotate("", xy=(5.6, 2.4), xytext=(3.6, 2.4), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6, mutation_scale=16))
    ax.text(4.6, 2.85, "stamp L\n(nodes 2, 3)", ha="center", fontsize=8.5)
    # the matrix
    x0, y0, cw, chh = 6.6, 0.4, 2.2, 1.2
    for i in range(3):
        ax.text(x0 + cw * (i + 0.5), y0 + 3 * chh + 0.35, str(i + 1), ha="center", fontsize=10, fontweight="bold")
        ax.text(x0 - 0.35, y0 + chh * (2.5 - i), str(i + 1), ha="center", va="center", fontsize=10, fontweight="bold")
        for j in range(3):
            hit = i >= 1 and j >= 1
            ax.add_patch(Rectangle((x0 + cw * j, y0 + chh * (2 - i)), cw, chh, fc="#FBE3D1" if hit else "white", ec="#999", lw=0.8))
            if hit:
                ax.text(x0 + cw * (j + 0.5), y0 + chh * (2.5 - i), ("+" if i == j else "−") + "%g" % gl, ha="center",
                        va="center", fontsize=10.5, color=ACCENT)
    save(fig, "stamp.svg")


def fig_ghost(v):
    fig, ax = canvas(2.6, 2.5, (-1.8, 3.4), (-1.3, 3.6))
    th = np.linspace(0, math.pi, 40)
    top = [(math.cos(t) * 1.0, 1.0 + math.sin(t) * 1.0) for t in th]
    bottom = []
    xs = np.linspace(-1.0, 1.0, 50)
    for x in xs:
        bottom.append((x, -0.6 + 0.18 * math.sin((x + 1) / 2 * 4 * 2 * math.pi)))
    pts = [(1.0, 1.0)] + top[1:] + [(-1.0, -0.6)] + bottom + [(1.0, 1.0)]
    ax.add_patch(Polygon(pts, closed=True, fc="white", ec=INK, lw=1.6))
    eyes(ax, 0, 1.3, 1.3)
    ax.add_patch(Ellipse((0, 0.5), 0.35, 0.45, fc=INK))
    bubble(ax, 1.6, 2.9, 3.0, 1.1, "±%.0f V...\nwoooo!" % abs(v), (0.8, 1.9), size=12)
    save(fig, "ghost.svg")


# ============================================================ the book
def box(kind, title, body):
    return "<div class='box %s'><div class='tag'>%s</div>%s</div>" % (kind, title, body)


def words(s):
    return "<p class='words'>In words: %s</p>" % s


def build():
    os.makedirs(FIG, exist_ok=True)
    print("figures ...")
    fig_cover(); fig_photo_movie(); fig_rl(); fig_trapezoid(); fig_companion(); fig_rl_companion()
    fig_network(); fig_rl_open(); fig_cast(); fig_pete_emma(); fig_loop()

    print("numbers ...")
    # ---------------- chapter 4: RL circuit
    dt = 1e-4
    rl = rl_circuit(dt).run(5e-3)
    check("rl.net", rl, ["R", "L"])
    exact = lambda t: 10 * (1 - np.exp(-t / 1e-3))
    rows = []
    for k, s in enumerate(rl[:10]):
        t = s["t"]
        rows.append((k + 1, t * 1e3, s["hist"]["L"], s["v"][2], s["i"]["L"], exact(t), exact(t - dt / 2)))
    s1, s2, s3 = rl[0], rl[1], rl[2]
    fig_sticky(s1["i"]["L"], s1["v"][2])

    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.6), gridspec_kw={"width_ratios": [1.5, 1]})
    tt = np.linspace(0, 5e-3, 500)
    ts = np.array([0] + [s["t"] for s in rl])
    isim = np.array([0] + [s["i"]["L"] for s in rl])
    for ax, xmax in ((axs[0], 5), (axs[1], 0.62)):
        ax.plot(tt * 1e3, exact(tt), color=GREY, lw=2.2, label="exact")
        ax.plot(ts * 1e3, isim, "o-", color=ACCENT, ms=3.2 if xmax < 1 else 2, lw=1.1, label="EMT, Δt = 0.1 ms")
        ax.set_xlim(0, xmax); ax.set_xlabel("time (ms)")
    axs[1].plot(tt * 1e3, exact(tt - dt / 2) * (tt >= dt / 2), "--", color=INK, lw=1, label="exact, Δt/2 later")
    axs[0].set_ylabel("current i (A)"); axs[0].set_ylim(0, 10.5); axs[1].set_ylim(0, 5.5)
    axs[0].legend(frameon=False, fontsize=9, loc="lower right"); axs[1].legend(frameon=False, fontsize=8.5, loc="upper left")
    axs[1].set_title("the first six steps", fontsize=10)
    fig.tight_layout(); save(fig, "rl_result.svg")

    study = []
    for d in (5e-4, 1e-4, 2e-5):
        r = rl_circuit(d).run(1e-3)
        study.append((d * 1e3, r[-1]["i"]["L"], abs(r[-1]["i"]["L"] - exact(1e-3)), abs(r[-1]["i"]["L"] - exact(1e-3 - d / 2))))

    # ---------------- chapter 5: network
    dtn = 5e-5
    net = network(dtn).run(0.08)
    check("network.net", net, ["R1", "L", "C", "LOAD"])
    G = net[0]["G"]
    GUU = G[1:, 1:]
    GUK = G[1:, :1]
    n1, n2 = net[0], net[1]
    gl, gc = dtn / (2 * 0.01), 2 * 10e-6 / dtn
    fig_stamp(gl)
    w = 2 * math.pi * 50
    ZL, ZC = 1j * w * 0.01, 1 / (1j * w * 10e-6)
    Zp = 1 / (1 / ZC + 1 / 100)
    V3 = 100 * Zp / (1 + ZL + Zp)
    fig, ax = plt.subplots(figsize=(6.6, 2.5))
    tn = np.array([s["t"] for s in net]) * 1e3
    ax.plot(tn, [s["v"][3] for s in net], color=ACCENT, lw=1.1, label="node 3, EMT")
    tt = np.linspace(0, 0.08, 2000)
    ax.plot(tt * 1e3, abs(V3) * np.cos(w * tt + np.angle(V3)), "--", color=INK, lw=0.9, label="steady state (phasor)")
    ax.set_xlabel("time (ms)"); ax.set_ylabel("voltage (V)"); ax.set_xlim(0, 80)
    ax.set_ylim(-125, 205)
    ax.legend(frameon=False, fontsize=9, loc="upper right", ncol=2)
    fig.tight_layout(); save(fig, "network_result.svg")
    v3 = np.array([s["v"][3] for s in net])
    peak = np.abs(v3).max()
    f0 = 1 / (2 * math.pi * math.sqrt(0.01 * 10e-6))

    # ---------------- chapter 6: fault
    netf = network(dtn, fault=(0.03, 0.06)).run(0.08)
    check("network_fault.net", netf, ["R1", "L", "C", "LOAD", "FAULT"])
    Gf_on = network(dtn, fault=(0.03, 0.06)).G(0.035)[1:, 1:]
    Gf_off = network(dtn, fault=(0.03, 0.06)).G(0.01)[1:, 1:]
    fig, axs = plt.subplots(2, 1, figsize=(6.6, 3.9), sharex=True)
    tf = np.array([s["t"] for s in netf]) * 1e3
    for ax in axs:
        ax.axvspan(30, 60, color=FAULT, lw=0)
    axs[0].plot(tf, [s["v"][3] for s in netf], color=ACCENT, lw=1.0)
    axs[0].set_ylabel("node 3 (V)")
    iL = np.array([s["i"]["L"] for s in netf])
    axs[1].plot(tf, iL, color=INK, lw=1.0)
    axs[1].set_ylabel("current in L (A)"); axs[1].set_xlabel("time (ms)"); axs[1].set_xlim(0, 80)
    axs[0].text(45, axs[0].get_ylim()[1] * 0.82, "fault", ha="center", color="#B03030", fontsize=10)
    fig.tight_layout(); save(fig, "fault_result.svg")
    ifault_peak = np.abs(iL[(tf > 30) & (tf < 60)]).max()

    # ---------------- chapter 8: opening the switch
    ro = rl_circuit(dt, opens=3e-3).run(5e-3)
    check("rl_open.net", ro, ["SWITCH", "R", "L"])
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.5), gridspec_kw={"width_ratios": [1.4, 1]})
    to = np.array([0] + [s["t"] for s in ro]) * 1e3
    io = np.array([0] + [s["i"]["L"] for s in ro])
    vo = np.array([0] + [s["v"][3] for s in ro])
    axs[0].plot(to, io, "-", color=INK, lw=1.3)
    axs[0].axvline(3, color=GREY, ls="--", lw=0.8)
    axs[0].set_xlabel("time (ms)"); axs[0].set_ylabel("current in L (A)"); axs[0].set_xlim(0, 5)
    axs[0].text(3.1, 5, "switch\nopens", fontsize=8.5, color="#444")
    m_ = (to >= 2.55) & (to <= 3.85)
    axs[1].plot(to[m_], vo[m_], "o-", color=ACCENT, ms=3.5, lw=0.9)
    axs[1].set_ylabel("voltage across L (V)")
    axs[1].axvline(3, color=GREY, ls="--", lw=0.8)
    axs[1].set_xlabel("time (ms)"); axs[1].set_title("zoomed in: every step", fontsize=10)
    fig.tight_layout(); save(fig, "rl_open_result.svg")
    k_open = int(round(3e-3 / dt)) - 1
    i_before = ro[k_open - 1]["i"]["L"]
    i_after = max(abs(ro[k_open + j]["i"]["L"]) for j in range(1, 20))
    v_before = ro[k_open - 1]["v"][3]
    v_after = [ro[k_open + j]["v"][3] for j in range(5)]
    fig_ghost(v_after[-1])

    print("html ...")
    H = []
    a = H.append
    r4 = rows[3]

    # ---------------------------------------------------------------- cover
    a("""
<section class="cover">
  <div class="titlebox">Tiny book of EMT simulation</div>
  <div class="subtitle">Closing a switch</div>
  <div class="tile-wrap"><div class="tile"><img src="fig/cover.svg" alt=""></div>
  <div class="tile reflection"><img src="fig/cover.svg" alt=""></div></div>
  <div class="author">By Amila Pathirana</div>
</section>
""")

    # ---------------------------------------------------------------- 1
    a("""
<h2>1. Introduction</h2>
<p>Flip a light switch. The light comes on. Boring?</p>
<p>Not for the network behind the switch. For a few thousandths of a second, its voltages and currents
swing, overshoot and ring like a bell before they calm down. Close a big circuit breaker, or let
lightning hit a line, and the same thing happens, only much more violently. These short, wild moments
are called <b>electromagnetic transients</b>. <b>EMT simulation</b> is how we calculate them on a
computer, so we can check that cables, transformers and breakers will survive them.</p>
<p>You do not need to be a maths wizard for this. If you know Ohm's law, you can do EMT simulation.
There is one clever trick in the middle, called the <i>companion model</i>, and this book will walk you
through it with real numbers, one step at a time.</p>
<p>The best way to learn it is with a simple example. But first, let us meet the cast.</p>

<h2>2. Meet the cast</h2>
<p>Every circuit in this book is built from three characters (Figure 1).</p>
<figure><img src="fig/cast.svg" style="width:100%">
<figcaption>Figure 1: The cast. The resistor has no memory; the inductor and the capacitor do.</figcaption></figure>
<p><b>The resistor</b> is the easy one. Its current depends only on its voltage <i>right now</i>:
<i>v</i> = <i>R i</i>. It does not care what happened a moment ago.</p>
<p><b>The inductor</b> hates sudden changes in its current. Its voltage is <i>L</i> times how fast its
current changes. To change the current suddenly, you would need an enormous voltage.</p>
<p><b>The capacitor</b> is the mirror image. It hates sudden changes in its voltage. Its current is
<i>C</i> times how fast its voltage changes.</p>
""" + box("remember", "Remember", """<p>The inductor and the capacitor have <b>memory</b>: what they do now depends
on what they did a moment ago. That memory is the whole difficulty of EMT simulation, and, as you will
see, also the whole trick.</p>"""))

    # ---------------------------------------------------------------- 3
    a("""
<h2>3. A photo and a movie</h2>
<p>Meet Phasor Pete and EMT Emma (Figure 2).</p>
<figure><img src="fig/pete_emma.svg" style="width:94%">
<figcaption>Figure 2: Pete has a photo. Emma has a movie.</figcaption></figure>
<p>Pete uses <b>phasors</b>, like most power-system programs (load flow, short circuit). A phasor is
like one photo of a sine wave, with a note on the back: "100 V, 50 Hz". That is all you need, as long
as the wave keeps going the same way forever.</p>
<p>But right after a switch closes, nothing is a steady sine wave yet. A photo cannot show what
happens next. For that, you need Emma's movie: one picture every small time step, Δt, where each
picture is worked out from the one before it (Figure 3).</p>
<figure><img src="fig/photo_movie.svg" style="width:100%">
<figcaption>Figure 3: A phasor is one photo. An EMT simulation is a movie with one frame every Δt.</figcaption></figure>
<p>A small Δt gives a smooth movie, but more frames to calculate. Typical time steps are 20 to 100 μs
for switching studies, a few μs for power electronics, and less than a microsecond for lightning.</p>
""")

    # ---------------------------------------------------------------- 4
    hist3 = s3["hist"]["L"]
    a(f"""
<h2>4. Our first EMT simulation</h2>
<p>Here is our circuit (Figure 4): a 10 V battery, a switch, a 1 Ω resistor and a 1 mH inductor.
Everything is at rest. At <i>t</i> = 0, the switch closes. What does the current <i>i</i> do?</p>
<figure><img src="fig/rl_circuit.svg" style="width:66%">
<figcaption>Figure 4: Our first circuit. Node 1 is the battery's terminal; node 2 is the top of the
inductor.</figcaption></figure>
<p>We picked this circuit for a reason: it is simple enough to solve exactly with pen and paper, so we
can check our simulation against the right answer:</p>
{eq("i(t) =", fr("V", "R"), br("1 −", sp("e", "− t / τ"), o="(", c=")"), ",&#8195;", "τ =", fr("L", "R"), "= 1 'ms")}
{words("the current climbs from 0 towards <i>V</i>/<i>R</i> = 10 A. After one time constant τ = 1 ms it has made it 63 % of the way.")}
<p>Our goal: get the same answer with nothing but Ohm's law and a calculator. No differential
equations.</p>

<h3>a. Chop time into small steps</h3>
<p>We pick a time step Δt = 0.1 ms, and only calculate the circuit at <i>t</i> = Δt, 2Δt, 3Δt, and
so on. The resistor is happy with that: at every instant, its current is just its voltage divided by
<i>R</i>. The inductor is the problem. Its rule links "now" to "a moment ago":</p>
{eq("v = L", fr("di", "dt"))}

<h3>b. The inductor's memory</h3>
<p>Turn the inductor's rule around, and add up its voltage over one time step:</p>
{eq("i(t) = i(t − Δt) +", fr("1", "L"), integral("t − Δt", "t", "v dt"), label=1)}
{words("current now = current one step ago + (1/<i>L</i>) × the area under the voltage curve during the step.")}
<p>We do not know exactly how the voltage changed during the step. So we guess the simplest thing: a
straight line from its old value to its new value. The area under a straight line is a trapezoid
(Figure 5). This guess is called the <b>trapezoidal rule</b>.</p>
<figure><img src="fig/trapezoid.svg" style="width:50%">
<figcaption>Figure 5: The trapezoidal rule. The area under the curve is replaced by the orange
trapezoid. The small grey sliver is our error.</figcaption></figure>
{eq(integral("t − Δt", "t", "v dt"), "≈", fr("Δt", "2"), br("v(t) + v(t − Δt)"), label=2)}
{words("area ≈ width × the average of the two heights.")}
<p>Put (2) into (1), and sort everything into two piles: things that belong to <i>now</i>, and things
that belong to the <i>last step</i>:</p>
{eq("i(t) =", fr("Δt", "2L"), "v(t) +", br("i(t − Δt) +", fr("Δt", "2L"), "v(t − Δt)"), label=3)}
<p>The bracket only holds things from the last step, which we already know. It is the inductor's
memory: a little note it writes to itself at the end of every step (Figure 6). We call it the
<b>history current</b>, <i>I</i><sub>hist</sub>.</p>
<figure><img src="fig/sticky.svg" style="width:84%">
<figcaption>Figure 6: The inductor's note to itself after step 1. These two numbers are all it needs
to remember.</figcaption></figure>

<h3>c. The disguise: the companion model</h3>
<p>"Companion model" sounds complex, but it is not. Look at equation (3) again:</p>
<ul>
<li>The first part, (Δt/2L) × <i>v</i>(<i>t</i>), is a conductance times the voltage now. That is
exactly what a resistor does, with resistance <i>R<sub>L</sub></i> = 2L/Δt.</li>
<li>The second part, <i>I</i><sub>hist</sub>, is already known. A current that is already known is
just a current source.</li>
</ul>
{eq("i(t) =", fr("v(t)", R_L), "+", I_HIST, "(t),&#8195;&#8195;", R_L, "=", fr("2L", "Δt"), ",&#8195;&#8195;", I_HIST, "(t) = i(t − Δt) +", fr("v(t − Δt)", R_L), label=4)}
<p>So for one time step, the inductor puts on a disguise: <b>a resistor in parallel with a current
source</b>. The resistor never changes (as long as Δt stays the same). Only the current source is
updated, once per step. The capacitor gets the same disguise, with <i>R<sub>C</sub></i> = Δt/2C and a
minus sign in its history current (Figure 7).</p>
<figure><img src="fig/companion.svg" style="width:86%">
<figcaption>Figure 7: The companion models. For one time step, an inductor or a capacitor between
nodes <i>k</i> and <i>m</i> becomes a resistor and a current source.</figcaption></figure>
<p>For our inductor: {m(R_L, "=", fr("2 × 0.001", "0.0001"), "= 20 'Ω")}.</p>

<h3>d. Solving one time step</h3>
<p>Swap the inductor for its disguise. Now the circuit has only a battery, two resistors and a current
source (Figure 8). Ohm's law can handle that.</p>
<figure><img src="fig/rl_companion.svg" style="width:70%">
<figcaption>Figure 8: The circuit of Figure 4 for one time step, with the inductor in disguise.</figcaption></figure>
<p>Node 1 is easy: the battery holds it at <i>V</i> = 10 V. Node 2 is unknown; call its voltage
<i>v</i><sub>2</sub>. Whatever current flows in through <i>R</i> must flow out through the companion
model (Kirchhoff's current law):</p>
{eq(fr(cat("V −", sb("v", "2")), "R"), "=", fr(sb("v", "2"), R_L), "+", I_HIST, "&#8195;⟹&#8195;",
    sb("v", "2"), "=", fr(cat(fr("V", "R"), "−", I_HIST), cat(fr("1", "R"), "+", fr("1", R_L))), label=5)}
{words("current in through <i>R</i> = current out through the disguised inductor. Solve for <i>v</i><sub>2</sub>.")}
<p>Then equation (4) gives the inductor current.</p>
<p>Stop and look at who is who in (5). On the left, <i>v</i><sub>2</sub> is <b>now</b>: the one thing
we do not know. On the right, <i>V</i>, <i>R</i> and <i>R<sub>L</sub></i> never change, and
<i>I</i><sub>hist</sub> is made entirely from the <b>last step</b>: its current and its voltage,
from the note. So every time step is the same little cycle (Figure 9): read the note from the
last step, solve for the voltages now, work out the currents now, write the note for the next step.</p>
<figure><img src="fig/loop.svg" style="width:100%">
<figcaption>Figure 9: One time step. The last step's voltages and currents live only in the notes
(the history currents); the solve uses them on the right-hand side and gives this step's voltages;
the currents follow, and become the last step of the next step.</figcaption></figure>
<p>Here are the four quantities of one inductor, and when each one is known:</p>
<table class="narrow">
<tr><th>symbol</th><th>where it comes from</th><th>when</th></tr>
<tr><td><i>v</i>(<i>t</i>−Δt), <i>i</i>(<i>t</i>−Δt)</td><td>the result of the last step</td><td>known before the step starts</td></tr>
<tr><td><i>I</i><sub>hist</sub>(<i>t</i>)</td><td>= <i>i</i>(<i>t</i>−Δt) + <i>v</i>(<i>t</i>−Δt) / <i>R<sub>L</sub></i></td><td>computed first, from the two above</td></tr>
<tr><td><i>v</i>(<i>t</i>)</td><td>from <i>G</i> <i>v</i>(<i>t</i>) = <i>i</i><sub>hist</sub>(<i>t</i>)</td><td>the unknown: the solve gives it</td></tr>
<tr><td><i>i</i>(<i>t</i>)</td><td>= <i>v</i>(<i>t</i>) / <i>R<sub>L</sub></i> + <i>I</i><sub>hist</sub>(<i>t</i>)</td><td>last, from <i>v</i>(<i>t</i>) and the note</td></tr>
</table>
<p>Put the last two rows together and you are back at equation (3), the trapezoidal rule:
<i>i</i>(<i>t</i>) − <i>i</i>(<i>t</i>−Δt) = (Δt/2<i>L</i>) [<i>v</i>(<i>t</i>) + <i>v</i>(<i>t</i>−Δt)]. So <i>G</i> <i>v</i> = <i>i</i> is
nothing but that rule, written for the whole circuit at once, with everything from <i>t</i>−Δt moved to the
right-hand side. The current <i>i</i>(<i>t</i>) never appears in it; it is read off afterwards, and
then it is the <i>i</i>(<i>t</i>−Δt) of the next step.</p>
""" + box("remember", "Remember", """<p>In <i>G</i> <i>v</i> = <i>i</i>, the two sides live at different times.
<b>Left: now</b>, the voltages we are solving for. <b>Right: the last step</b>, folded into the
history currents. The currents of this step are not on the right-hand side; they are computed
afterwards from the new voltages, and only then do they go into the next right-hand side.</p>""") + f"""


<h3>e. Step by step, with numbers</h3>
<p>Before <i>t</i> = 0, everything is at rest: <i>i</i> = 0 and <i>v</i><sub>2</sub> = 0. Grab a
calculator and follow along. I tried not to skip a single step.</p>
<p class="step">Step 1, <i>t</i> = 0.1 ms</p>
<p>Read the note from the last step (everything was zero):</p>
{eq(I_HIST, "= 0 +", fr("0", "20"), "=", num(0, 4), "'A")}
{eq(sb("v", "2"), "=", fr("10 / 1 − 0", "1 / 1 + 1 / 20"), "=", fr("10", "1.05"), "=", num(s1['v'][2]), "'V")}
{eq("i =", fr(num(s1['v'][2]), "20"), "+ 0 =", num(s1['i']['L']), "'A")}
<p class="step">Step 2, <i>t</i> = 0.2 ms</p>
<p>Now the note says what happened in step 1:</p>
{eq(I_HIST, "=", num(s1['i']['L']), "+", fr(num(s1['v'][2]), "20"), "=", num(s2['hist']['L']), "'A")}
{eq(sb("v", "2"), "=", fr(cat("10 −", num(s2['hist']['L'])), "1.05"), "=", num(s2['v'][2]), "'V")}
{eq("i =", fr(num(s2['v'][2]), "20"), "+", num(s2['hist']['L']), "=", num(s2['i']['L']), "'A")}
<p class="step">Step 3, <i>t</i> = 0.3 ms</p>
{eq(I_HIST, "=", num(s2['i']['L']), "+", fr(num(s2['v'][2]), "20"), "=", num(hist3), "'A")}
{eq(sb("v", "2"), "=", fr(cat("10 −", num(hist3)), "1.05"), "=", num(s3['v'][2]), "'V")}
{eq("i =", fr(num(s3['v'][2]), "20"), "+", num(hist3), "=", num(s3['i']['L']), "'A")}
""" + box("recipe", "The recipe", """<p>Every step is the same three lines, over and over:</p>
<ol><li><b>Read the note:</b> work out the history current from the last step.</li>
<li><b>Solve:</b> find the node voltage with Ohm's law.</li>
<li><b>Write a new note:</b> work out the current, and keep it for the next step.</li></ol>
<p>That is all a real EMT program does, just a few million times faster than you.</p>""")
      + box("try", "Try it", f"""<p>Work out step 4 yourself, using the step 3 numbers. Then check your answer in
Table 1.</p><p class="answer">Answer: <i>I</i><sub>hist</sub> = {f(r4[2])} A, <i>v</i><sub>2</sub> =
{f(r4[3])} V, <i>i</i> = {f(r4[4])} A.</p>""") + """
<table>
<caption>Table 1: The first ten steps, next to the exact current from the formula.</caption>
<tr><th>step</th><th><i>t</i> (ms)</th><th><i>I</i><sub>hist</sub> (A)</th><th><i>v</i><sub>2</sub> (V)</th>
<th><i>i</i> (A)</th><th>exact <i>i</i> (A)</th></tr>
""" + "".join("<tr><td>%d</td><td>%.1f</td><td class='num'>%s</td><td class='num'>%s</td><td class='num'>%s</td><td>%s</td></tr>"
              % (r[0], r[1], f(r[2]), f(r[3]), f(r[4]), f(r[5])) for r in rows) + f"""
</table>

<h3>f. How good is the answer?</h3>
<p>Figure 10 puts our simulation next to the exact answer. After a few steps the two curves almost
overlap, but the simulation is always a little behind. Why?</p>
<figure><img src="fig/rl_result.svg" style="width:100%">
<figcaption>Figure 10: Simulated (orange dots) and exact (grey) current. Right: the first six steps.
The simulation matches the exact curve shifted by half a time step (black dashes).</figcaption></figure>
<p>Our simulation has slow reflexes. In step 1, the trapezoidal rule assumed the inductor voltage rose
in a straight line from 0 V (when everything was at rest) to {f(s1['v'][2], 2)} V. In reality it jumped
to 10 V the instant the switch closed. A voltage that ramps up over one step acts, on average, like a
voltage that jumps halfway through the step. So the simulation behaves as if the switch had closed half
a step late, at Δt/2. Shift the exact curve by Δt/2 and the two agree very well (right of Figure 10).
A smaller Δt means smaller reflex delay:</p>
<table>
<caption>Table 2: The current at <i>t</i> = 1 ms for three time steps (exact answer: {f(exact(1e-3))} A).</caption>
<tr><th>Δt (ms)</th><th>simulated <i>i</i> (A)</th><th>error (A)</th><th>error after the Δt/2 shift (A)</th></tr>
""" + "".join("<tr><td>%g</td><td class='num'>%s</td><td>%s</td><td>%s</td></tr>"
              % (r[0], f(r[1]), f(r[2], 5), f(r[3], 5)) for r in study) + """
</table>
""" + box("remember", "Remember", """<p><b>Δt must be small compared with the fastest thing you want to see.</b>
Here the circuit's time constant is 1 ms, so a step of 0.1 ms is plenty.</p>"""))

    # ---------------------------------------------------------------- 5
    a(f"""
<h2>5. A bigger network</h2>
<p>Writing Kirchhoff's law by hand for every node gets boring fast. Computers do it with a matrix.
Figure 11 shows a slightly bigger circuit: a 100 V, 50 Hz source feeding a line (<i>R</i>1 and
<i>L</i>), a cable that adds some capacitance <i>C</i>, and a 100 Ω load. The time step is
Δt = 50 μs.</p>
<figure><img src="fig/network.svg" style="width:80%">
<figcaption>Figure 11: The network of chapters 5 and 6. Nodes 1, 2 and 3; node 0 is ground.</figcaption></figure>

<h3>a. Everything becomes a conductance</h3>
<p>In disguise, every element is a conductance (1 / resistance), plus maybe a current source:</p>
<table class="narrow">
<tr><th>element</th><th>between nodes</th><th>conductance</th></tr>
<tr><td><i>R</i>1 = 1 Ω</td><td>1 and 2</td><td class='num'>1/1 = 1 S</td></tr>
<tr><td><i>L</i> = 10 mH</td><td>2 and 3</td><td class='num'>Δt/2L = {gl:g} S</td></tr>
<tr><td><i>C</i> = 10 μF</td><td>3 and ground</td><td class='num'>2C/Δt = {gc:g} S</td></tr>
<tr><td>LOAD = 100 Ω</td><td>3 and ground</td><td class='num'>1/100 = 0.01 S</td></tr>
</table>

<h3>b. The rubber stamp</h3>
<p>Now we build the conductance matrix <b>G</b>, one element at a time, with a rubber stamp
(Figure 12). The stamp for a conductance <i>g</i> between nodes <i>k</i> and <i>m</i> adds +<i>g</i> to
the two diagonal spots (<i>k</i>,<i>k</i>) and (<i>m</i>,<i>m</i>), and −<i>g</i> to the two
off-diagonal spots (<i>k</i>,<i>m</i>) and (<i>m</i>,<i>k</i>). If one end is ground, only the
diagonal spot of the other end gets stamped.</p>
<figure><img src="fig/stamp.svg" style="width:88%">
<figcaption>Figure 12: Stamping the inductor (between nodes 2 and 3) into <b>G</b>.</figcaption></figure>
<p>Stamp all four elements, and you get:</p>
{eq(BG, "=", mat(G, 4))}
<p>Kirchhoff's current law for all the nodes at once is then one neat line:</p>
{eq(BG, BV, "(t) =", sb(BI, "hist"), "(t)", label=6)}
{words("conductances × the voltages <b>now</b> = the history currents written at the <b>last step</b>.")}
<p><b>v</b>(<i>t</i>) is this step's unknown. <b>i</b><sub>hist</sub>(<i>t</i>) is built before the solve,
from the last step's voltages and currents: a history current <i>I</i><sub>hist</sub> flowing from node
<i>k</i> to node <i>m</i> takes <i>I</i><sub>hist</sub> out of entry <i>k</i> and puts it into entry <i>m</i>.
The branch currents of this step are not in it; they come afterwards, element by element, from
<i>i</i> = <i>v</i>/<i>R</i><sub>eq</sub> + <i>I</i><sub>hist</sub>, and feed the next step's right-hand side.
(In the rest of the chapter we write just <b>i</b> for this vector.)</p>

<h3>c. The boss node</h3>
<p>The source is the boss of node 1: it tells node 1 exactly what voltage to have. So node 1 is
<b>known</b> (K), and nodes 2 and 3 are <b>unknown</b> (U). Split <b>G</b> into those parts:</p>
{eq(sb(BG, "UU"), "=", mat(GUU, 4), ",&#8195;&#8195;", sb(BG, "UK"), "=", mat(GUK, 0))}
<p>and solve only for the unknown voltages:</p>
{eq(sb(BG, "UU"), sb(BV, "U"), "=", sb(BI, "U"), "−", sb(BG, "UK"), sb(BV, "K"), label=7)}
{words("move everything we already know (the boss's voltage) to the right-hand side, then solve for the rest.")}
<p>Here is the beautiful part: <b>G</b> stays the same from one step to the next. Only the right-hand
side changes. So the computer does the hard work (factorising <b>G</b><sub>UU</sub>, with an LU or
LDL<sup>T</sup> factorisation in EMTSim) only once. After that, each step is quick. That is why an EMT
program can run millions of steps.</p>

<h3>d. One step, with numbers</h3>
<p class="step">Step 1, <i>t</i> = 50 μs</p>
<p>The boss sets node 1 to {m("100 cos(2π × 50 × 0.00005) =", num(n1['vK'][0]), "'V")}. Everything was
at rest before, so every note is blank: all history currents are zero, and <b>i</b><sub>U</sub> = 0.</p>
{eq(sb(BI, "U"), "−", sb(BG, "UK"), sb(BV, "K"), "=", col([0, 0], 0), "−", mat(GUK, 0),
    mat([[n1['vK'][0]]], 4), "=", col(n1['rhs'], 4))}
{eq(sb(BV, "U"), "=", sp(sb(BG, "UU"), "−1"), col(n1['rhs'], 4), "=", col(n1['v'][2:], 4))}
<p>Node 2 jumps to {f(n1['v'][2], 2)} V, while node 3 barely moves: the inductor is not letting much
current through yet. Its current is {m("0.0025 × (", num(n1['v'][2]), "−", num(n1['v'][3]), ") =",
                                      num(n1['i']['L'], 6), "'A")}, and the capacitor's is
{m(num(0.4, 1), "×", num(n1['v'][3]), "=", num(n1['i']['C'], 6), "'A")}.</p>
<p class="step">Step 2, <i>t</i> = 100 μs</p>
<p>First, read the notes (equation 4 and Figure 7):</p>
{eq(sb("I", "hist,L"), "=", num(n1['i']['L'], 6), "+", num(gl, 4), "×", num(n1['v'][2] - n1['v'][3]), "=", num(n2['hist']['L'], 6), "'A")}
{eq(sb("I", "hist,C"), "= −", num(n1['i']['C'], 6), "−", num(gc, 1), "×", num(n1['v'][3]), "=", num(n2['hist']['C'], 6), "'A")}
<p><i>I</i><sub>hist,L</sub> flows from node 2 to node 3, and <i>I</i><sub>hist,C</sub> from node 3 to
ground. Then solve:</p>
{eq(sb(BI, "U"), "=", col([n2['inj'][1], n2['inj'][2]], 6), ",&#8195;", sb(BV, "K"), "=", num(n2['vK'][0]),
    "'V,&#8195;", sb(BV, "U"), "=", col(n2['v'][2:], 4))}
<p>And so on, step after step. Figure 13 shows the movie of node 3 for the first 80 ms. When the
source switches on, the inductor and the capacitor play catch with the energy, at their natural
frequency {m("1 / (2π", "<msqrt>" + cat("LC") + "</msqrt>", ") ≈", num(f0, 0), "'Hz")}. The voltage
overshoots to {f(peak, 1)} V. Then the load calms the ringing down, and we end up exactly where Phasor
Pete's photo says we should.</p>
<figure><img src="fig/network_result.svg" style="width:100%">
<figcaption>Figure 13: Voltage of node 3 after the source switches on, and Pete's steady state
(a phasor calculation of the same circuit).</figcaption></figure>
""" + box("recipe", "The recipe of every EMT program", """
<p><b>Once:</b> put every element in disguise, stamp <b>G</b>, and factorise <b>G</b><sub>UU</sub>.</p>
<p><b>Every step</b> (<i>t</i> = <i>t</i> + Δt):</p>
<ol>
<li>Did a switch or a fault change? Then re-stamp <b>G</b> and factorise it again.</li>
<li>Ask the bosses (the sources) for the known voltages <b>v</b><sub>K</sub>.</li>
<li>Read every inductor's and capacitor's note: the history currents. Add them up into <b>i</b>.</li>
<li>Solve <b>G</b><sub>UU</sub> <b>v</b><sub>U</sub> = <b>i</b><sub>U</sub> − <b>G</b><sub>UK</sub> <b>v</b><sub>K</sub>.</li>
<li>Work out every element's current, write the new notes, and save the results.</li>
</ol>
<p>This is exactly the order EMTSim uses.</p>"""))

    # ---------------------------------------------------------------- 6
    a(f"""
<h2>6. Switches and faults</h2>
<p>A switch is just a resistor with two moods: tiny when it is closed (EMTSim uses 0.01 Ω) and huge when
it is open (1 MΩ). A fault is a switch from a node to ground.</p>
<p>When a switch changes its mood, its conductance changes, so <b>G</b> changes. That is step 1 of the
recipe: re-stamp <b>G</b>, factorise it again, and carry on. The notes (history currents) are not
touched, so the inductors and capacitors remember their currents and voltages across the switching
instant, just like real ones.</p>
<p>Let us short-circuit node 3 to ground at 30 ms and remove the short at 60 ms (Figure 14).</p>
<figure><img src="fig/network_fault.svg" style="width:80%">
<figcaption>Figure 14: The network of Figure 11, with a fault from node 3 to ground.</figcaption></figure>
<p>While the fault is open, it stamps a tiny 1/10<sup>6</sup> S onto <b>G</b><sub>33</sub>. When it
closes, it stamps 1/0.01 = 100 S:</p>
{eq("before:&#8195;", sb(BG, "UU"), "=", mat(Gf_off, 6), "&#8195;&#8195;during:&#8195;", sb(BG, "UU"), "=", mat(Gf_on, 4))}
<p>Figure 15 shows what happens. When the fault hits, node 3 collapses to almost zero, and the capacitor
dumps its charge into the fault. The current in <i>L</i> grows to about {f(ifault_peak, 1)} A; the
inductor is the only thing limiting it. When the fault is removed, the capacitor and the inductor play
catch again at about {f0:.0f} Hz, until the network settles back to normal.</p>
<figure><img src="fig/fault_result.svg" style="width:100%">
<figcaption>Figure 15: Voltage of node 3 and current in the inductor, with the fault (shaded) from 30 to
60 ms.</figcaption></figure>
""")

    # ---------------------------------------------------------------- 7
    def netfile(name):
        with open(os.path.join(NETS, name), encoding="utf-8") as fh:
            return fh.read().strip().replace("&", "&amp;").replace("<", "&lt;")
    a(f"""
<h2>7. Try it yourself with EMTSim</h2>
<p>Every number in this book was calculated with the same recipe EMTSim uses, and checked against it.
EMTSim gives the same numbers, to the last digit. So you can run every example yourself.</p>
<p>An EMTSim netlist has one line per element: its type, the two nodes it connects (0 is ground) and its
value. This is the circuit of chapter 4:</p>
<pre>{netfile("rl.net")}</pre>
<p>and this is the network with the fault, from chapter 6:</p>
<pre>{netfile("network_fault.net")}</pre>
<p>Run them from the EMTSim folder, with either the C++ or the Python solver:</p>
<pre>build\\emtsim.exe docs\\book\\circuits\\rl.net
python main.py docs\\book\\circuits\\network_fault.net
python plot.py</pre>
<p>The node voltages go to <code>results.csv</code>, and the element currents to
<code>currents.csv</code>. Compare the first lines of <code>results.csv</code> with Table 1. Or open
the circuit editor with <code>gui.bat</code>, draw the circuit, press Run, and drag across a plot to
zoom in on the first few steps.</p>
""" + box("try", "Try it", """<ul style="margin:0">
<li>Make Δt ten times bigger. How big does the error in Table 2 get?</li>
<li>Move the fault to node 2. Is the fault current bigger or smaller? Why?</li>
<li>Change the load to 10 Ω. Does the ringing die out faster or slower?</li></ul>"""))

    # ---------------------------------------------------------------- 8
    a(f"""
<h2>8. The ghost in the machine</h2>
<p>So far we have only closed switches. Let us take the circuit of chapter 4 again, wait until the
current has almost reached 10 A, and then <i>open</i> the switch at <i>t</i> = 3 ms (Figure 16).</p>
<figure><img src="fig/rl_open_circuit.svg" style="width:60%">
<figcaption>Figure 16: The circuit of chapter 4, with a switch that opens at 3 ms.</figcaption></figure>
<p>Remember the inductor's personality: it hates sudden changes in its current. In a real circuit, when
the switch opens, the inductor fights back by pushing the voltage across the opening contacts up very
high (that is why breakers draw arcs), and the current dies away very quickly. Now look at what our
simulation says (Figure 17).</p>
<figure><img src="fig/rl_open_result.svg" style="width:100%">
<figcaption>Figure 17: Left: current in the inductor. Right: voltage across the inductor at every step
around the opening. Both calculated with the trapezoidal rule and Δt = 0.1 ms.</figcaption></figure>
<img src="fig/ghost.svg" class="ghost" alt="">
<p>The current looks fine. Just before the switch opens, it is {f(i_before, 2)} A, and from the next
step on, it is practically zero (less than {i_after:.4f} A).</p>
<p>But look at the voltage across the inductor. Just before the opening, it is {f(v_before, 2)} V. In
the steps after it, the simulation says {", ".join(f(x, 1) for x in v_after)} V, and so on. It flips
sign at every single step, for as long as we care to simulate. No current flows, and yet the inductor
voltage keeps jumping between plus and minus {f(abs(v_after[-1]), 0)} V. Spooky.</p>
<p>Does the voltage in a real circuit do that? Or is there a ghost in our simulation? Think about it
before you turn the page; chapter 9 catches it.</p>
""" + box("try", "Ghost hunting", """<p>Use equation (4), with the 1 MΩ of the open switch in series with
<i>R<sub>L</sub></i> = 20 Ω. What happens to <i>I</i><sub>hist</sub> from one step to the next? Who
wrote the wrong note?</p>"""))

    # ---------------------------------------------------------------- 9: chatter removal and interpolation
    cda, half = rl_open_cda(dt)
    check("rl_open.net", cda, ["R", "L"], "--chatter=on", "--interp=on")
    k3 = int(round(3e-3 / dt)) - 1            # index of the 3.0 ms step
    s30, s31, s32 = cda[k3], cda[k3 + 1], cda[k3 + 2]
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.5), gridspec_kw={"width_ratios": [1.4, 1]})
    tc = np.array([0] + [x["t"] for x in cda]) * 1e3
    vc = np.array([0] + [x["v"][3] for x in cda])
    axs[0].plot(to, vo, ".-", color=GREY, lw=0.8, ms=3, label="plain (chapter 8)")
    axs[0].plot(tc, vc, "o-", color=ACCENT, ms=3, lw=1.2, label="with chatter removal")
    axs[0].axvline(3, color=GREY, ls="--", lw=0.8)
    axs[0].set_xlabel("time (ms)"); axs[0].set_ylabel("voltage across L (V)"); axs[0].set_xlim(0, 5)
    axs[0].legend(frameon=False, fontsize=8.5, loc="lower left")
    m_ = (tc >= 2.55) & (tc <= 3.65)
    axs[1].plot(tc[m_], vc[m_], "o-", color=ACCENT, ms=4, lw=1.2)
    axs[1].plot([half["t"] * 1e3], [half["v"][3]], "o", mfc="white", mec=ACCENT, ms=6)
    axs[1].annotate("half step\n(not recorded)", (half["t"] * 1e3, half["v"][3]), (3.25, half["v"][3] * 0.55), fontsize=8,
                    arrowprops=dict(arrowstyle="-", color="#666", lw=0.7))
    axs[1].axvline(3, color=GREY, ls="--", lw=0.8)
    axs[1].set_xlabel("time (ms)"); axs[1].set_title("around the opening", fontsize=10)
    fig.tight_layout(); save(fig, "cda_result.svg")

    # interpolation: the demo circuit, run with EMTSim (dt = 0.2 ms, switch closes at 1.02 ms)
    demo = os.path.join(ROOT, "circuits", "demo_interpolation.net")
    late = emtsim(demo)                      # as in the file: interpolation off, chatter removal on
    exact_ip = lambda t: 10 * (1 - np.exp(-(t - 1.02e-3) / 1e-3)) * (t > 1.02e-3)
    have_demo = late is not None
    if have_demo:
        onv, oni = emtsim(demo, "--interp=on")
        tl, il = late[0][:, 0], late[1]["L"]
        ti, ii = onv[:, 0], oni["L"]
        fig, ax = plt.subplots(figsize=(6.6, 2.5))
        tt = np.linspace(0, 6e-3, 600)
        ax.plot(tt * 1e3, exact_ip(tt), color=GREY, lw=2.2, label="exact")
        ax.plot(np.r_[0, tl] * 1e3, np.r_[0, il], "s--", color=INK, ms=3.5, lw=0.9, label="switch at the next time point (1.2 ms)")
        ax.plot(np.r_[0, ti] * 1e3, np.r_[0, ii], "o-", color=ACCENT, ms=3.5, lw=1.1, label="interpolation: switch at 1.02 ms")
        ax.set_xlim(0.8, 3.2); ax.set_ylim(-0.3, 9); ax.set_xlabel("time (ms)"); ax.set_ylabel("current i (A)")
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
        fig.tight_layout(); save(fig, "interp_result.svg")
        k14 = int(np.argmin(np.abs(tl - 1.4e-3)))
        i_late, i_on, i_ex = il[k14], ii[k14], exact_ip(1.4e-3)
    hist_half = s30['i']['L']                 # the half step's history current: the current alone
    a(f"""
<h2>9. Catching the ghost</h2>
<p>Look at the inductor's note again, equation (4): <i>I</i><sub>hist</sub> = <i>i</i>(<i>t</i>−Δt) +
<i>v</i>(<i>t</i>−Δt)/<i>R<sub>L</sub></i>. It remembers the <b>voltage</b> of the last step as well as the
current. After the switch opens the current must be zero at every step, so equation (4) says
0 = <i>v</i>(<i>t</i>)/<i>R<sub>L</sub></i> + <i>v</i>(<i>t</i>−Δt)/<i>R<sub>L</sub></i>, that is</p>
{eq("v(t) = − v(t − Δt)")}
<p>The voltage has to flip sign every step, forever. That is the ghost: the trapezoidal rule averages
the voltage over each step, so a sudden jump can only be represented by an endless alternation around
the true value. Nothing in the circuit does this; the note does. Programs call it <b>numerical
chatter</b>, and it appears after every switching that jumps an inductor voltage or a capacitor
current.</p>

<h3>a. Chatter removal: write a shorter note</h3>
<p>The cure is to forget the old voltage, just for a moment. Right after a switching we take two
<b>half steps</b> of Δt/2 with a different rule, <i>backward Euler</i>, whose note holds only the
current:</p>
{eq("i(t) = i(t − Δt/2) +", fr("Δt", "2L"), "v(t),&#8195;&#8195;", I_HIST, "= i(t − Δt/2)")}
{words("the same conductance Δt/2L as before, so <b>G</b> does not change; only the note is shorter.")}
<p>With a step of exactly Δt/2, backward Euler has the same conductance as the trapezoidal rule, so
nothing has to be rebuilt or refactorised. After the two half steps the trapezoidal rule takes over
again as if nothing had happened. This is called <i>critical damping adjustment</i>, CDA; the
half-step solution is not written to the output, so the time grid stays the same.</p>
<p class="step">The circuit of chapter 8, with chatter removal</p>
<p>At 3.0 ms the switch is still closed in that sample: <i>i</i> = {f(s30['i']['L'])} A,
<i>v</i><sub>3</sub> = {f(s30['v'][3])} V. Then it opens. First half step, to 3.05 ms: the note is the current
alone, {f(hist_half)} A. The open switch (1 MΩ) lets almost nothing through, so the inductor can only
keep its note by pushing its voltage to</p>
{eq(sb("v", "3"), "(3.05 'ms) ≈ −", R_L, "×", num(hist_half), "=", num(half['v'][3], 1), "'V,&#8195;", "i ≈", num(half['i']['L'], 4), "'A")}
<p>That is the inductor's one protest at being switched off (in a real breaker it is the arc). Second
half step, to 3.1 ms: the note is now {f(half['i']['L'], 4)} A, so <i>v</i><sub>3</sub> = {f(s31['v'][3], 3)} V and
<i>i</i> = {f(s31['i']['L'], 4)} A. From here the trapezoidal rule continues with a note of
{f(s31['i']['L'] + s31['v'][3] / 20, 4)} A: at 3.2 ms <i>v</i><sub>3</sub> = {f(s32['v'][3], 3)} V. A few millivolts remain, the echo of the
{f(half['i']['L'], 4)} A that leaked through the 1 MΩ during the half step; the 187 V ghost is gone (Figure 18).</p>
<figure><img src="fig/cda_result.svg" style="width:100%">
<figcaption>Figure 18: The inductor voltage of chapter 8 with and without chatter removal. The hollow
point is the half-step solution, which is not part of the output.</figcaption></figure>

<h3>b. Interpolation: switch at the right instant</h3>
<p>The second problem with a fixed time step is simpler to see. A switch is told to close at 1.02 ms,
but the simulation only visits 1.0 ms and 1.2 ms. If it can only act at a time point, it closes at
1.2 ms, and the whole current curve starts almost a step late (Figure 19, squares).</p>
<p>The fix: take the step from 1.0 to 1.2 ms as usual, then notice that the switch should have acted
at a fraction α = (1.02 − 1.0) / 0.2 = 0.1 of the way through it. Go back to that instant by
<b>interpolating the whole state</b> (every node voltage, every note) between the two time points,</p>
{eq("state(1.02 'ms) = (1 − α) state(1.0 'ms) + α state(1.2 'ms)")}
<p>close the switch there, rebuild <b>G</b>, and continue from 1.02 ms: two half steps (chatter
removal, since the switching jumps the inductor voltage), then normal steps. The simulation is now on a
time grid shifted by 0.02 ms, so after the next step it interpolates back onto the original grid;
the output has the same time points as before. The same trick finds a breaker's current zero, or the
instant a diode's current crosses zero, inside a step.</p>
""" + (f"""
<figure><img src="fig/interp_result.svg" style="width:100%">
<figcaption>Figure 19: A 10 V battery switched onto 1 Ω + 1 mH at 1.02 ms, with Δt = 0.2 ms (the
circuit <code>demo_interpolation.net</code>). At 1.4 ms the exact current is {f(i_ex, 2)} A; the switch
acting at the next time point gives {f(i_late, 2)} A, interpolation {f(i_on, 2)} A.</figcaption></figure>
""" if have_demo else "") + """
<p>The two belong together. Interpolation alone would still leave the trapezoidal rule's half-step
delay after the switching (and the chatter); chatter removal alone acts one step late. With both, a
switching acts at its instant as accurately as the trapezoidal rule itself allows.</p>
""" + box("try", "Try it", """<p>Both methods are on by default in EMTSim; the book's circuits switch them off with
the line <code>.options chatter=off interp=off</code>. Run chapter 8's circuit with
<code>run.bat circuits\\demo_chatter.net --chatter=on</code>, or open <code>demo_chatter</code> and
<code>demo_interpolation</code> in the editor: Run, tick <b>Hold run</b>, tick the option under Run
settings, Run again, and the two results are drawn over each other. <code>docs/switching.md</code> has the
details and where the code is.</p>"""))

    # ---------------------------------------------------------------- cheat sheet
    a(f"""
<h2 class="newpage">Cheat sheet</h2>
<div class="cheat">
<div><h4>The cast</h4>
{eq("v = R i,&#8195;&#8195;v = L", fr("di", "dt"), ",&#8195;&#8195;i = C", fr("dv", "dt"))}</div>
<div><h4>Inductor in disguise (trapezoidal rule)</h4>
{eq("i(t) =", fr("v(t)", R_L), "+", I_HIST, "(t),&#8195;", R_L, "=", fr("2L", "Δt"))}
{eq(I_HIST, "(t) = i(t − Δt) +", fr("v(t − Δt)", R_L))}</div>
<div><h4>Capacitor in disguise</h4>
{eq("i(t) =", fr("v(t)", R_C), "+", I_HIST, "(t),&#8195;", R_C, "=", fr("Δt", "2C"))}
{eq(I_HIST, "(t) = − i(t − Δt) −", fr("v(t − Δt)", R_C))}</div>
<div><h4>The stamp of a conductance <i>g</i> between nodes <i>k</i> and <i>m</i></h4>
{eq(sb(BG, "kk"), "<mo>+=</mo>", "g,&#8195;", sb(BG, "mm"), "<mo>+=</mo>", "g,&#8195;", sb(BG, "km"), "<mo>−=</mo>", "g,&#8195;", sb(BG, "mk"), "<mo>−=</mo>", "g")}</div>
<div><h4>One time step (left: the voltages now; right: history currents from the last step)</h4>
{eq(sb(BG, "UU"), sb(BV, "U"), "=", sb(BI, "U"), "−", sb(BG, "UK"), sb(BV, "K"))}</div>
<div><h4>Choosing Δt</h4><p style="text-indent:0;text-align:center">much smaller than the fastest thing you
want to see: the simulation acts about Δt/2 late after a sudden jump.</p></div>
<div><h4>Switches and faults</h4><p style="text-indent:0;text-align:center">closed 0.01 Ω, open 1 MΩ.
When one changes: re-stamp <b>G</b>, factorise again, keep the notes.</p></div>
<div><h4>Chatter removal (after a switching: two backward-Euler half steps, same <b>G</b>)</h4>
{eq(I_HIST, "= i(t − Δt/2)", "&#8195;(L),&#8195;&#8195;", I_HIST, "= −", fr("v(t − Δt/2)", R_C), "&#8195;(C)")}</div>
<div><h4>Interpolation (a switching at fraction α of the step)</h4>
{eq("state(t + αΔt) = (1 − α) state(t) + α state(t + Δt)")}</div>
</div>
""")

    html = ("<!doctype html><html lang='en'><head><meta charset='utf-8'><title>Tiny book of EMT simulation</title>"
            "<style>%s</style></head><body>%s</body></html>" % (CSS, "\n".join(H)))
    out = os.path.join(HERE, "emt_book.html")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(html)
    print("  wrote", os.path.relpath(out, ROOT))
    return out


CSS = r"""
@page { size: letter; margin: 0.95in 1.05in 0.9in 1.05in;
  @bottom-left { content: counter(page) "  |  Page"; font: 10pt Calibri, 'Segoe UI', sans-serif; color: #7f7f7f;
                 border-top: 0.6pt solid #d0d0d0; vertical-align: top; padding-top: 4pt; width: 100%; } }
@page :first { @bottom-left { content: none; border: none; } }
html { background: white; }
body { font-family: Georgia, 'Times New Roman', serif; font-size: 12.5pt; line-height: 1.32; color: #111;
       margin: 0; background: white; }
math { font-family: 'Cambria Math', 'STIX Two Math', serif; font-size: 1.05em; }
.num, mn.num { color: #C55A11; }
p { margin: 0 0 0.75em; text-align: justify; text-indent: 1.6em; hyphens: auto; }
ul, ol { margin: 0.2em 0 0.8em; }
li { margin-bottom: 0.25em; }
h2 { background: #000; color: #fff; font-size: 15.5pt; font-weight: bold; padding: 2pt 10pt 1pt;
     margin: 1.3em 0 0.6em; break-after: avoid; }
h2.newpage { break-before: page; margin-top: 0; }
h3 { background: #000; color: #fff; font-size: 14pt; font-weight: bold; padding: 1pt 10pt 0;
     margin: 1.2em 0 0.55em 0.9in; break-after: avoid; }
h2 + p, h3 + p { break-before: avoid; }
p.step { text-indent: 0; font-weight: bold; display: inline-block; background: #d9d9d9; padding: 0 5pt;
         margin: 0.6em 0 0.4em 0.9in; }
p.words { text-indent: 0; text-align: center; font-style: italic; color: #444; font-size: 11.5pt; margin: -0.2em 0 0.9em; }
.eq { position: relative; margin: 0.5em 0 0.7em; break-inside: avoid; text-align: center; }
.eqn { position: absolute; right: 0; top: 50%; transform: translateY(-50%); }
figure { margin: 0.9em 0 1em; text-align: center; break-inside: avoid; }
figure img { display: block; margin: 0 auto; }
figcaption { text-align: left; font-size: 11.5pt; margin-top: 0.35em; }
table { border-collapse: collapse; margin: 0.6em auto 1em; font-size: 11pt; break-inside: avoid; }
caption { caption-side: top; text-align: left; font-size: 11.5pt; padding-bottom: 0.4em; }
th, td { border-bottom: 0.6pt solid #bbb; padding: 2pt 9pt; text-align: right; }
th { border-bottom: 1.2pt solid #000; font-weight: bold; }
td:first-child, th:first-child { text-align: left; }
pre { font-family: Consolas, 'Courier New', monospace; font-size: 9.6pt; background: #f2f2f2;
      border-left: 3pt solid #000; padding: 6pt 10pt; white-space: pre-wrap; break-inside: avoid; line-height: 1.25; }
code { font-family: Consolas, monospace; font-size: 0.9em; }
/* boxes */
.box { border: 1.3pt solid #000; border-radius: 7pt; padding: 0 12pt 4pt; margin: 1.1em 0 1.1em; break-inside: avoid; }
.box .tag { display: inline-block; background: #000; color: #fff; font-family: 'Segoe Print', cursive; font-size: 11pt;
            padding: 0 9pt 1pt; border-radius: 0 0 6pt 6pt; margin-bottom: 5pt; }
.box p { text-indent: 0; }
.box.remember { border-color: #C55A11; }
.box.remember .tag { background: #C55A11; }
.box.try { background: #f4f4f4; }
.box.recipe { background: #FFF8D6; }
p.answer { font-size: 10.5pt; color: #555; border-top: 0.6pt dashed #999; padding-top: 3pt; }
img.ghost { float: right; width: 2.5in; margin: 0 0 0.3em 0.8em; }
/* cheat sheet */
.cheat > div { border: 1pt solid #000; border-radius: 6pt; padding: 2pt 10pt 2pt; margin-bottom: 8pt; break-inside: avoid; }
.cheat h4 { margin: 3pt 0 0; font-size: 11pt; }
.cheat .eq { margin: 0.3em 0 0.4em; }
/* cover */
.cover { height: 9.1in; display: flex; flex-direction: column; align-items: center; break-after: page; }
.titlebox { background: #000; color: #fff; font-family: 'Segoe Script', 'Segoe Print', cursive; font-size: 40pt;
            line-height: 1.25; text-align: center; padding: 18pt 28pt 22pt; width: 82%; margin-top: 0.1in; }
.subtitle { font-family: 'Segoe Script', 'Segoe Print', cursive; color: #7f7f7f; font-size: 30pt; margin: 0.35in 0 0.25in; }
.tile-wrap { display: flex; flex-direction: column; align-items: center; }
.tile { width: 1.65in; height: 3.75in; border-radius: 10pt; background: white; border: 1pt solid #cfcfcf;
        box-shadow: inset 0 0 14pt rgba(0,0,0,0.25), 0 1pt 3pt rgba(0,0,0,0.2); overflow: hidden;
        display: flex; align-items: center; justify-content: center; }
.tile img { width: 88%; height: 92%; object-fit: contain; }
.tile.reflection { height: 0.9in; transform: scaleY(-1); opacity: 0.35; margin-top: 3pt; border: none; box-shadow: none;
                   -webkit-mask-image: linear-gradient(to bottom, transparent 0%, rgba(0,0,0,0.9) 100%);
                   mask-image: linear-gradient(to bottom, transparent 0%, rgba(0,0,0,0.9) 100%);
                   align-items: flex-end; }
.tile.reflection img { height: 3.45in; width: 88%; object-fit: contain; object-position: bottom; }
.author { margin-top: auto; background: #000; color: #fff; font-family: 'Segoe Script', 'Segoe Print', cursive;
          font-size: 21pt; width: 100%; text-align: center; padding: 3pt 0 5pt; }
@media screen { body { max-width: 6.4in; margin: 0.6in auto; } .cover { height: auto; min-height: 9.1in; } }
"""


def print_pdf(html):
    candidates = [r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
                  r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
                  r"C:\Program Files\Google\Chrome\Application\chrome.exe",
                  shutil.which("msedge") or "", shutil.which("google-chrome") or "", shutil.which("chromium") or ""]
    browser = next((c for c in candidates if c and os.path.exists(c)), None)
    if not browser:
        print("  no Edge or Chrome found: open emt_book.html and print it to PDF")
        return
    pdf = os.path.join(HERE, "TinyBookOfEMT.pdf")
    prof = tempfile.mkdtemp()
    try:
        subprocess.run([browser, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                        "--user-data-dir=" + prof, "--virtual-time-budget=5000",
                        "--print-to-pdf=" + pdf, "file:///" + html.replace("\\", "/")],
                       check=True, capture_output=True, timeout=120)
    finally:
        shutil.rmtree(prof, ignore_errors=True)
    print("  wrote", os.path.relpath(pdf, ROOT))


if __name__ == "__main__":
    out = build()
    if "--no-pdf" not in sys.argv:
        print("pdf ...")
        print_pdf(out)
