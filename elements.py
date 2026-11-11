"""Circuit elements of the Python EMTSim (same models as src/elements.cpp).

Every element connects k terminals (node numbers, 0 = ground). For each time
step it behaves, as seen from the network, like

    i = g v + h

where v are the terminal voltages, i the currents flowing from each terminal
node into the element, g a k x k conductance matrix (stamp) and h a vector of
history currents (history).

Per time step the network calls, in this order:
    time_event(t)   scheduled state changes (switches)        -> True: g changed
    history(t)      history currents for the step
    ... network solution ...
    check(t, v)     state changes that depend on the solution  -> True: g changed,
                    (diodes, arresters, saturation)                 re-solve the step
    update(t, v)    commit the solution

Ideal voltage sources (and machines' internal EMFs) fix the voltage of some
terminals: known_terminals() / known_voltages(t); the network then tells them
the current they deliver with set_known_currents().

Chatter removal and switching interpolation add: a backward-Euler half step when config.euler
is set (same conductances, no voltage term in the histories; config.step is the step length),
save_state/load_state for the continuous state, and for switching elements begin_step /
event_alpha (where in the step, 0..1, it should have switched; None: not) / apply_event.
"""
import bisect
import cmath
import functools
import math

import numpy as np

import config

PHASE = ("a", "b", "c")


def _phase_name(name, k, n):
    return name + "." + PHASE[k] if n == 3 else name if n == 1 else name + "." + str(k + 1)


def _get(obj, name):
    return functools.reduce(getattr, name.split("."), obj)


def _set(obj, name, value):
    *path, last = name.split(".")
    setattr(functools.reduce(getattr, path, obj), last, value)


def _frac(x):
    return min(1.0, max(0.0, x))


class Element:
    nonlinear = False              # True if check() can change the state
    has_events = False             # True for switches, faults, diodes and thyristors
    valve = False                  # diodes and thyristors: switched by events, not by check()
    state_vars = ()                # attributes holding the continuous state (numbers, arrays, complex)

    def __init__(self, type_, nodes):
        self.type = type_
        self.name = ""
        self.nodes = list(nodes)
        self.idx = None            # matrix index of each terminal (-1 = ground), set by the network

    def size(self):
        return len(self.nodes)

    def stamp(self):               # k x k conductance
        raise NotImplementedError

    def history(self, t):          # k history currents
        raise NotImplementedError

    def time_event(self, t):
        return False

    def check(self, t, v):
        return False

    def update(self, t, v):
        raise NotImplementedError

    def known_terminals(self):     # positions in nodes[] whose voltage the element fixes
        return []

    def known_voltages(self, t):
        return []

    def set_known_currents(self, i):
        pass

    def current_names(self):
        return [self.name]

    def currents(self):
        raise NotImplementedError

    def save_state(self, out):
        for name in self.state_vars:
            v = _get(self, name)
            if isinstance(v, np.ndarray):
                out.extend(v.ravel().tolist())
            elif isinstance(v, complex):
                out.extend((v.real, v.imag))
            else:
                out.append(float(v))

    def load_state(self, vals, pos):
        for name in self.state_vars:
            v = _get(self, name)
            if isinstance(v, np.ndarray):
                n = v.size
                _set(self, name, np.array(vals[pos:pos + n], dtype=float).reshape(v.shape))
                pos += n
            elif isinstance(v, complex):
                _set(self, name, complex(vals[pos], vals[pos + 1]))
                pos += 2
            else:
                _set(self, name, float(vals[pos]))
                pos += 1
        return pos

    def begin_step(self, a, b):
        pass

    def event_alpha(self, a, dt):  # fraction 0..1 of the step a -> a+dt where the element switches, or None
        return None

    def apply_event(self, t):
        return False


# ---------------------------------------------------------------- two-terminal branches

class TwoTerminal(Element):
    """i (from strnode to stpnode) = V / Reff + ihistory,  V = v[0] - v[1]"""

    def __init__(self, type_, strnode, stpnode):
        super().__init__(type_, [strnode, stpnode])
        self.Reff = 1.0
        self.Ilast = 0.0
        self.Vlast = 0.0
        self.ihistory = 0.0

    @property
    def strnode(self):
        return self.nodes[0]

    @property
    def stpnode(self):
        return self.nodes[1]

    @property
    def brnType(self):
        return self.type

    def stamp(self):
        g = 1.0/self.Reff
        return np.array([[g, -g], [-g, g]])

    def calc_history(self, t):
        self.ihistory = 0.0

    def history(self, t):
        self.calc_history(t)
        return [self.ihistory, -self.ihistory]

    def update(self, t, v):
        self.Vlast = v[0] - v[1]
        self.Ilast = self.Vlast/self.Reff + self.ihistory

    def currents(self):
        return [self.Ilast]

    def save_state(self, out):
        out.append(self.Ilast)
        out.append(self.Vlast)

    def load_state(self, vals, pos):
        self.Ilast, self.Vlast = vals[pos], vals[pos + 1]
        return pos + 2


class RLC(TwoTerminal):
    """Resistor, inductor or capacitor with trapezoidal-rule companion model"""

    def __init__(self, brnType=None, strnode=None, stpnode=None, value=None):
        super().__init__(brnType, strnode, stpnode)
        if not value > 0:
            raise ValueError("%s: value must be positive" % brnType)
        if brnType == "R":
            self.Reff = value
        elif brnType == "L":
            self.Reff = 2.0*value/config.Dt         # 2L/Dt
        elif brnType == "C":
            self.Reff = config.Dt/(2.0*value)       # Dt/2C
        else:
            raise ValueError("RLC: type must be R, L or C")

    def calc_history(self, t):
        if self.type == "L":
            self.ihistory = self.Ilast if config.euler else self.Ilast + self.Vlast/self.Reff
        elif self.type == "C":
            self.ihistory = -self.Vlast/self.Reff if config.euler else -self.Ilast - self.Vlast/self.Reff
        else:
            self.ihistory = 0.0


class Source(TwoTerminal):
    """Sinusoidal voltage source E(t) = magnitude cos(2 pi f t + angle), raising stpnode above strnode.
    resistance > 0: Norton equivalent;  resistance None or 0: ideal (one terminal on ground)"""

    def __init__(self, brnType="S", strnode=None, stpnode=None, magnitude=None, angle=None, frequency=None,
                 resistance=None):
        super().__init__("V", strnode, stpnode)
        self.magnitude = magnitude
        self.angle = angle
        self.frequency = frequency
        self.ideal = not resistance or resistance <= 0
        if self.ideal and strnode != 0 and stpnode != 0:
            raise ValueError("An ideal voltage source must have one node connected to ground (node 0)")
        self.Reff = None if self.ideal else resistance
        self.E = 0.0

    def Sourceupdate(self, TheTime):
        self.E = self.e(TheTime)

    def e(self, t):
        return self.magnitude*math.cos(2.0*math.pi*self.frequency*t + self.angle*math.pi/180.0)

    def stamp(self):
        return np.zeros((2, 2)) if self.ideal else super().stamp()

    def history(self, t):
        if self.ideal:
            return [0.0, 0.0]
        self.ihistory = self.e(t)/self.Reff      # Norton equivalent: E/R in parallel with R
        return [self.ihistory, -self.ihistory]

    def update(self, t, v):
        if self.ideal:
            self.Vlast = v[0] - v[1]
        else:
            super().update(t, v)

    def known_terminals(self):
        if not self.ideal:
            return []
        return [1 if self.nodes[0] == 0 else 0]

    def known_voltages(self, t):
        return [self.e(t) if self.nodes[0] == 0 else -self.e(t)]

    def set_known_currents(self, i):
        self.Ilast = i[0] if self.nodes[0] == 0 else -i[0]


class Switch(TwoTerminal):
    """Time-controlled switch (breaker): Ron closed, Roff open; toggles at each time in t_switch.
    With opens_at_current_zero, an opening command takes effect at the next current zero
    (or as soon as |i| <= chop, to model current chopping)."""

    def __init__(self, brnType="SW", strnode=None, stpnode=None, Ron=0.01, Roff=1.0e6, closed=True, t_switch=(),
                 opens_at_current_zero=False, chop=0.0):
        super().__init__(brnType, strnode, stpnode)
        self.Ron = Ron
        self.Roff = Roff
        self.closed0 = closed
        self.t_switch = sorted(t_switch)
        self.closed = closed
        self.Reff = Ron if closed else Roff
        self.opens_at_current_zero = opens_at_current_zero
        self.chop = chop
        self._pending_open = False
        self._zero_seen = False
        self._t_cmd = 0.0
        self._ia = 0.0

    has_events = True

    @staticmethod
    def _tol():                # time accumulates round-off, so switch times are matched with a tolerance
        return 1e-6*config.Dt

    def state_at(self, TheTime):
        n = sum(1 for ts in self.t_switch if TheTime >= ts - self._tol())
        return self.closed0 != (n % 2 == 1)

    def describe(self, closed):
        return "closed" if closed else "opened"

    def label(self):
        return "switch"

    def _set_state(self, closed, t):   # switch, and forget any pending open command
        self.closed = closed
        self.Reff = self.Ron if closed else self.Roff
        self._pending_open = self._zero_seen = False
        print("t = %.5f s: %s %s (nodes %d-%d)" % (t, self.label(), self.describe(closed), self.nodes[0], self.nodes[1]))

    def _command_open(self, t):        # wait for the next current zero
        self._pending_open, self._zero_seen, self._t_cmd = True, False, t
        print("t = %.5f s: %s open command (nodes %d-%d), waiting for current zero"
              % (t, self.label(), self.nodes[0], self.nodes[1]))

    def time_event(self, t):
        want = self.state_at(t)
        if want == self.closed:
            if want:
                self._pending_open = False
            return False
        if want or not self.opens_at_current_zero:   # closing is immediate; so is opening without current-zero detection
            self._set_state(want, t)
            return True
        if not self._pending_open:
            self._command_open(t)
        if not self._zero_seen:
            return False
        self._set_state(False, t)
        return True

    def update(self, t, v):
        iprev = self.Ilast
        super().update(t, v)
        if self._pending_open and self.closed and (iprev*self.Ilast <= 0.0 or abs(self.Ilast) <= self.chop):
            self._zero_seen = True

    # ---- events (chatter removal / interpolation)
    def begin_step(self, a, b):
        self._ia = self.Ilast
        if not self.opens_at_current_zero:
            return
        want = self.state_at(b)
        if self.closed and not want and not self._pending_open:
            self._command_open(max([a] + [ts for ts in self.t_switch if ts <= b + self._tol()]))
        elif want and self._pending_open:
            self._pending_open = False

    def event_alpha(self, a, dt):
        tol = self._tol()
        best = None
        for ts in self.t_switch:                     # scheduled closing / opening
            if a + tol < ts <= a + dt + tol:
                want = self.state_at(ts)
                if want != self.closed and (want or not self.opens_at_current_zero):
                    best = _frac((ts - a)/dt)
                    break
        if self._pending_open and self.closed:       # breaker: current zero (or chopping level)
            ia, ib = self._ia, self.Ilast
            al = None
            if ia == 0.0:
                al = 0.0
            elif ia*ib <= 0.0:
                al = ia/(ia - ib)
            if self.chop > 0.0 and abs(ib) <= self.chop:
                ac = 0.0 if abs(ia) <= self.chop else (abs(ia) - self.chop)/(abs(ia) - abs(ib))
                al = ac if al is None else min(al, ac)
            if al is not None:
                al = _frac(max(al, (self._t_cmd - a)/dt))
                best = al if best is None else min(best, al)
        return best

    def apply_event(self, t):
        want = self.state_at(t)
        if want == self.closed or (not want and self.opens_at_current_zero and not self._pending_open):
            return False
        self._set_state(want, t)
        return True


class Fault(Switch):
    """Fault: closed from t_apply until t_clear (t_clear None = permanent)"""

    def __init__(self, brnType="F", strnode=None, stpnode=0, Ron=0.01, Roff=1.0e6, t_apply=None, t_clear=None):
        times = [t_apply] if t_clear is None else [t_apply, t_clear]
        super().__init__(brnType=brnType, strnode=strnode, stpnode=stpnode, Ron=Ron, Roff=Roff, closed=False, t_switch=times)
        self.t_apply = t_apply
        self.t_clear = t_clear

    def describe(self, closed):
        return "applied" if closed else "cleared"

    def label(self):
        return "fault"


# ---------------------------------------------------------------- coupled R-L

def sequence_matrix(x1, x0):
    """Symmetric components -> 3x3 phase matrix with self (x0 + 2 x1)/3 and mutual (x0 - x1)/3"""
    s, m = (x0 + 2.0*x1)/3.0, (x0 - x1)/3.0
    return np.array([[s, m, m], [m, s, m], [m, m, s]])


class RLBlock:
    """n coupled series branches: i = A v + A [v_prev + (2L/dt - R) i_prev],  A = (R + 2L/dt)^-1"""

    def __init__(self, R, L):
        R = np.atleast_2d(np.asarray(R, dtype=float))
        L = np.atleast_2d(np.asarray(L, dtype=float))
        self.n = R.shape[0]
        self.A = np.linalg.inv(R + (2.0/config.Dt)*L)
        self.B = (2.0/config.Dt)*L - R
        self.L2 = (2.0/config.Dt)*L
        self.vprev = np.zeros(self.n)
        self.iprev = np.zeros(self.n)
        self.ih = np.zeros(self.n)

    def stamp(self, g, frm, to):
        for i in range(self.n):
            for j in range(self.n):
                a = self.A[i, j]
                g[frm[i], frm[j]] += a
                g[to[i], to[j]] += a
                g[frm[i], to[j]] -= a
                g[to[i], frm[j]] -= a

    def history(self, h, frm, to):
        self.ih = self.A @ (self.L2 @ self.iprev) if config.euler else self.A @ (self.vprev + self.B @ self.iprev)
        for i in range(self.n):
            h[frm[i]] += self.ih[i]
            h[to[i]] -= self.ih[i]

    def branch_voltage(self, v, frm, to):
        return np.array([v[frm[i]] - v[to[i]] for i in range(self.n)])

    def update(self, vb):
        self.iprev = self.A @ vb + self.ih
        self.vprev = vb


class CoupledRL(Element):
    """n coupled series R-L branches from[k] -> to[k] (nodes = from, then to)"""

    state_vars = ("rl.vprev", "rl.iprev")

    def __init__(self, frm, to, R, L):
        super().__init__("RL", list(frm) + list(to))
        n = len(frm)
        self._frm, self._to = list(range(n)), list(range(n, 2*n))
        self.rl = RLBlock(R, L)

    @classmethod
    def from_sequence(cls, frm, to, R1, L1, R0, L0):
        return cls(frm, to, sequence_matrix(R1, R0), sequence_matrix(L1, L0))

    def stamp(self):
        g = np.zeros((self.size(), self.size()))
        self.rl.stamp(g, self._frm, self._to)
        return g

    def history(self, t):
        h = np.zeros(self.size())
        self.rl.history(h, self._frm, self._to)
        return h

    def update(self, t, v):
        self.rl.update(self.rl.branch_voltage(v, self._frm, self._to))

    def current_names(self):
        return [_phase_name(self.name, k, self.rl.n) for k in range(self.rl.n)]

    def currents(self):
        return list(self.rl.iprev)


class CoupledC(Element):
    """n coupled shunt capacitances to ground (capacitance matrix C)"""

    state_vars = ("vprev", "iprev")

    def __init__(self, nodes, C):
        super().__init__("C", nodes)
        self.G = (2.0/config.Dt)*np.atleast_2d(np.asarray(C, dtype=float))
        n = len(nodes)
        self.vprev = np.zeros(n)
        self.iprev = np.zeros(n)
        self.ih = np.zeros(n)

    @classmethod
    def from_sequence(cls, nodes, C1, C0):
        return cls(nodes, sequence_matrix(C1, C0))

    def stamp(self):
        return self.G.copy()

    def history(self, t):
        self.ih = -(self.G @ self.vprev) if config.euler else -(self.G @ self.vprev + self.iprev)
        return self.ih

    def update(self, t, v):
        v = np.asarray(v, dtype=float)
        self.iprev = self.G @ v + self.ih
        self.vprev = v

    def current_names(self):
        return [_phase_name(self.name, k, self.size()) for k in range(self.size())]

    def currents(self):
        return list(self.iprev)


# ---------------------------------------------------------------- Bergeron line

class BergeronLine(Element):
    """Travelling-wave line with M modes (1: single-phase; 3: transposed three-phase line, Clarke
    modal transformation). Per mode: surge impedance Zc, travel time tau, total series resistance R
    (lumped R/4, R/2, R/4 as in the EMTP Bergeron model). nodes = k end (M), then m end (M)."""

    def __init__(self, k_end, m_end, Zc, tau, R):
        super().__init__("TL", list(k_end) + list(m_end))
        M = len(k_end)
        if len(m_end) != M or len(Zc) != M or len(tau) != M or len(R) != M or M not in (1, 3):
            raise ValueError("Bergeron line: 1 or 3 phases, with Zc, tau and R per mode")
        self.M = M
        if M == 1:
            self.T = np.array([[1.0]])
        else:   # orthonormal Clarke transformation: mode 0 (zero sequence), modes 1 and 2 (aerial)
            a, b, c, d = 1/math.sqrt(3), math.sqrt(2/3), 1/math.sqrt(6), 1/math.sqrt(2)
            self.T = np.array([[a, a, a], [b, -c, -c], [0.0, d, -d]])
        self.tau = list(tau)
        self.Z, self.h = [], []
        for j in range(M):
            if not Zc[j] > 0:
                raise ValueError("Bergeron line: surge impedance must be positive")
            if tau[j] < config.Dt*(1 - 1e-9):
                raise ValueError("Bergeron line: travel time must be at least one time step")
            self.Z.append(Zc[j] + R[j]/4.0)
            self.h.append((Zc[j] - R[j]/4.0)/(Zc[j] + R[j]/4.0))
        self.Yp = self.T.T @ np.diag([1.0/z for z in self.Z]) @ self.T
        # modal histories, one sample per solution point at times ts (sample 0 = start time)
        self.ts = [config.srtTime]
        self.vk = [[0.0] for _ in range(M)]
        self.ik = [[0.0] for _ in range(M)]
        self.vm = [[0.0] for _ in range(M)]
        self.im = [[0.0] for _ in range(M)]
        self.Ihk = np.zeros(M)
        self.Ihm = np.zeros(M)
        self.ik_phase = np.zeros(M)
        self.im_phase = np.zeros(M)

    def stamp(self):
        M = self.M
        g = np.zeros((2*M, 2*M))
        g[:M, :M] = self.Yp
        g[M:, M:] = self.Yp
        return g

    def _delayed(self, s, t, tau):
        x = t - tau
        ts = self.ts
        if x <= ts[0]:
            return 0.0 if x < ts[0] else s[0]
        i = bisect.bisect_right(ts, x) - 1
        if i + 1 >= len(ts):
            return s[-1]
        f = (x - ts[i])/(ts[i + 1] - ts[i])
        return s[i]*(1.0 - f) + s[i + 1]*f

    def history(self, t):
        for j in range(self.M):
            Z, hh, tau = self.Z[j], self.h[j], self.tau[j]
            vkd, ikd = self._delayed(self.vk[j], t, tau), self._delayed(self.ik[j], t, tau)
            vmd, imd = self._delayed(self.vm[j], t, tau), self._delayed(self.im[j], t, tau)
            self.Ihk[j] = 0.5*(1 + hh)*(-vmd/Z - hh*imd) + 0.5*(1 - hh)*(-vkd/Z - hh*ikd)
            self.Ihm[j] = 0.5*(1 + hh)*(-vkd/Z - hh*ikd) + 0.5*(1 - hh)*(-vmd/Z - hh*imd)
        return np.concatenate([self.T.T @ self.Ihk, self.T.T @ self.Ihm])

    def update(self, t, v):
        M = self.M
        vkm = self.T @ np.asarray(v[:M], dtype=float)
        vmm = self.T @ np.asarray(v[M:], dtype=float)
        ikm = np.empty(M)
        imm = np.empty(M)
        for j in range(M):
            ikm[j] = vkm[j]/self.Z[j] + self.Ihk[j]
            imm[j] = vmm[j]/self.Z[j] + self.Ihm[j]
            self.vk[j].append(vkm[j]); self.ik[j].append(ikm[j])
            self.vm[j].append(vmm[j]); self.im[j].append(imm[j])
        self.ts.append(t)
        self.ik_phase = self.T.T @ ikm
        self.im_phase = self.T.T @ imm

    # state: number of samples, time and values of the newest sample, phase currents. Loading a
    # state between two saved ones (a fractional count) replaces the newest sample.
    def save_state(self, out):
        out.append(float(len(self.ts)))
        out.append(self.ts[-1])
        for j in range(self.M):
            out.extend((self.vk[j][-1], self.ik[j][-1], self.vm[j][-1], self.im[j][-1]))
        out.extend(self.ik_phase.tolist())
        out.extend(self.im_phase.tolist())

    def load_state(self, vals, pos):
        # keep k samples, then the stored sample becomes the newest: appended if the count is fractional
        # (a state interpolated between two time points), otherwise it replaces sample k-1
        c = vals[pos]
        k = int(math.floor(c + 1e-9))
        keep = k if c - k > 1e-9 else k - 1
        for x in [self.ts] + self.vk + self.ik + self.vm + self.im:
            del x[keep:]
        self.ts.append(vals[pos + 1])
        for j in range(self.M):
            q = pos + 2 + 4*j
            self.vk[j].append(vals[q]); self.ik[j].append(vals[q + 1])
            self.vm[j].append(vals[q + 2]); self.im[j].append(vals[q + 3])
        pos += 2 + 4*self.M
        self.ik_phase = np.array(vals[pos:pos + self.M], dtype=float)
        self.im_phase = np.array(vals[pos + self.M:pos + 2*self.M], dtype=float)
        return pos + 2*self.M

    def current_names(self):
        return ([_phase_name(self.name + ".k", p, self.M) for p in range(self.M)] +
                [_phase_name(self.name + ".m", p, self.M) for p in range(self.M)])

    def currents(self):
        return list(self.ik_phase) + list(self.im_phase)


# ---------------------------------------------------------------- transformer

class TransformerCore(Element):
    """Ideal transformer (ratio n = V1/V2) with series leakage R, L on the primary side.
    nodes = p, q (primary), s, t (secondary). Branch voltage vb = (vp - vq) - n (vs - vt).
    The magnetising branch is a separate element across p-q."""

    state_vars = ("vprev", "i")

    def __init__(self, p, q, s, t, n, R, L):
        super().__init__("XF", [p, q, s, t])
        if not n > 0:
            raise ValueError("Transformer: ratio must be positive")
        if not R + 2.0*L/config.Dt > 0:
            raise ValueError("Transformer: leakage R and L cannot both be zero")
        self.a = np.array([1.0, -1.0, -n, n])
        self.n = n
        self.A = 1.0/(R + 2.0*L/config.Dt)
        self.B = 2.0*L/config.Dt - R
        self.L2 = 2.0*L/config.Dt
        self.vprev = 0.0
        self.i = 0.0
        self.ih = 0.0

    def stamp(self):
        return self.A*np.outer(self.a, self.a)

    def history(self, t):
        self.ih = self.A*self.L2*self.i if config.euler else self.A*(self.vprev + self.B*self.i)
        return self.a*self.ih

    def update(self, t, v):
        vb = float(self.a @ np.asarray(v, dtype=float))
        self.i = self.A*vb + self.ih
        self.vprev = vb

    def current_names(self):
        return [self.name + ".1", self.name + ".2"]

    def currents(self):
        return [self.i, self.n*self.i]


# ---------------------------------------------------------------- non-linear elements

def _stamp2(g):
    return np.array([[g, -g], [-g, g]])


class SaturableInductor(Element):
    """Non-linear inductor with a symmetric piecewise-linear flux-current curve through the origin:
    points (i_k, lambda_k), increasing; beyond the last point the last slope continues."""

    nonlinear = True
    state_vars = ("flux", "vprev", "i")

    def __init__(self, frm, to, i_points, flux_points):
        super().__init__("LSAT", [frm, to])
        if not i_points or len(i_points) != len(flux_points):
            raise ValueError("LSAT: give pairs of current and flux")
        self.ip, self.lp = [0.0], [0.0]
        for i, l in zip(i_points, flux_points):
            if not (i > self.ip[-1] and l > self.lp[-1]):
                raise ValueError("LSAT: currents and fluxes must be positive and increasing")
            self.ip.append(i)
            self.lp.append(l)
        self.seg, self.sign = 0, 1
        self.flux = 0.0
        self.vprev = 0.0
        self.i = 0.0
        self.ih = 0.0

    def _slope(self, seg):
        k = min(seg, len(self.ip) - 2)
        return (self.lp[k + 1] - self.lp[k])/(self.ip[k + 1] - self.ip[k])

    def _segment_of(self, lam):
        a = abs(lam)
        sign = 1 if lam >= 0 else -1
        last = len(self.lp) - 1
        for k in range(last):
            if a <= self.lp[k + 1]:
                return k, sign
        return last, sign

    def stamp(self):
        return _stamp2(config.Dt/(2.0*self._slope(self.seg)))

    def _flux0(self):
        # flux before this step's own voltage: lambda_prev + dt/2 v_prev (trapezoidal), lambda_prev (Euler)
        return self.flux if config.euler else self.flux + 0.5*config.Dt*self.vprev

    def history(self, t):
        # segment: i = s i0 + (lambda - s lambda0) / L, with lambda = flux0 + dt/2 v
        L, s = self._slope(self.seg), self.sign
        self.ih = s*self.ip[self.seg] + (self._flux0() - s*self.lp[self.seg])/L
        return [self.ih, -self.ih]

    def check(self, t, v):
        lam = self._flux0() + 0.5*config.Dt*(v[0] - v[1])
        seg, sign = self._segment_of(lam)
        if seg != self.seg or (seg != 0 and sign != self.sign):
            self.seg, self.sign = seg, sign
            return True
        return False

    def update(self, t, v):
        vb = v[0] - v[1]
        self.flux = self._flux0() + 0.5*config.Dt*vb
        self.i = config.Dt/(2.0*self._slope(self.seg))*vb + self.ih
        self.vprev = vb

    def currents(self):
        return [self.i]


class SurgeArrester(Element):
    """Metal-oxide surge arrester: i = Iref (|v| / Vref)^q sign(v), as a piecewise-linear resistor"""

    nonlinear = True
    state_vars = ("i",)

    def __init__(self, frm, to, Vref, Iref=1000.0, q=26.0):
        super().__init__("MOV", [frm, to])
        if not (Vref > 0 and Iref > 0 and q >= 1):
            raise ValueError("MOV: need Vref > 0, Iref > 0, q >= 1")
        # break points from 0.5 Vref to 2.5 Vref, the current growing by 30 % per segment
        self.vp, self.ip = [0.0], [0.0]
        ratio = 1.3**(1.0/q)
        x = 0.5
        while x < 2.5*ratio:
            self.vp.append(x*Vref)
            self.ip.append(Iref*x**q)
            x *= ratio
        self.seg, self.sign = 0, 1
        self.i = 0.0

    def _g(self, seg):
        k = min(seg, len(self.vp) - 2)
        return (self.ip[k + 1] - self.ip[k])/(self.vp[k + 1] - self.vp[k])

    def _i0(self, seg, sign):
        return sign*(self.ip[seg] - self._g(seg)*self.vp[seg])

    def _segment_of(self, v):
        a = abs(v)
        sign = 1 if v >= 0 else -1
        last = len(self.vp) - 1
        for k in range(last):
            if a <= self.vp[k + 1]:
                return k, sign
        return last, sign

    def stamp(self):
        return _stamp2(self._g(self.seg))

    def history(self, t):
        c = self._i0(self.seg, self.sign)
        return [c, -c]

    def check(self, t, v):
        seg, sign = self._segment_of(v[0] - v[1])
        if seg != self.seg or (seg != 0 and sign != self.sign):
            self.seg, self.sign = seg, sign
            return True
        return False

    def update(self, t, v):
        self.i = self._g(self.seg)*(v[0] - v[1]) + self._i0(self.seg, self.sign)

    def currents(self):
        return [self.i]


class Diode(Element):
    """Diode (anode, cathode): Ron (in series with Vf) when on, Roff when off.
    Turns on when v > Vf, off when the current reverses."""

    nonlinear = True
    has_events = True
    valve = True
    state_vars = ("i", "v")

    def __init__(self, anode, cathode, Ron=1e-3, Roff=1e6, Vf=0.0):
        super().__init__("D", [anode, cathode])
        if not (Ron > 0 and Roff > Ron):
            raise ValueError("Diode: need 0 < Ron < Roff")
        self.Ron, self.Roff, self.Vf = Ron, Roff, Vf
        self.on = False
        self.i = 0.0
        self.v = 0.0
        self._ia = self._va = 0.0

    def gate(self, t):
        return True

    def stamp(self):
        return _stamp2(1.0/self.Ron if self.on else 1.0/self.Roff)

    def history(self, t):
        c = -self.Vf/self.Ron if self.on else 0.0       # on: i = (v - Vf) / Ron
        return [c, -c]

    def check(self, t, v):
        vb = v[0] - v[1]
        if not self.on:
            if vb > self.Vf and self.gate(t):
                self.on = True
                return True
            return False
        if (vb - self.Vf)/self.Ron < 0.0:                # current reversed
            self.on = False
            return True
        return False

    def update(self, t, v):
        vb = v[0] - v[1]
        self.i = (vb - self.Vf)/self.Ron if self.on else vb/self.Roff
        self.v = vb

    def currents(self):
        return [self.i]

    # ---- events: turn on where the voltage crosses Vf (with the gate fired), off where the current crosses zero
    def begin_step(self, a, b):
        self._ia, self._va = self.i, self.v

    def event_alpha(self, a, dt):
        b = a + dt
        if not self.on:
            vb = self.v
            if vb <= self.Vf:
                return None
            va = self._va
            al = 0.0 if va >= self.Vf else (self.Vf - va)/(vb - va)
            t = a + al*dt
            if not self.gate(t):
                if not self.gate(b):
                    return None
                lo, hi = t, b                   # the gate fires during the step: find when
                for _ in range(50):
                    mid = 0.5*(lo + hi)
                    if self.gate(mid):
                        hi = mid
                    else:
                        lo = mid
                al = (hi - a)/dt
            return _frac(al)
        ib = self.i
        if ib >= 0.0:
            return None
        ia = self._ia
        if ia <= 0.0:
            return 0.0
        return _frac(ia/(ia - ib))

    def apply_event(self, t):
        self.on = not self.on
        return True


class Thyristor(Diode):
    """Thyristor: a diode that can only turn on while its gate is fired. Simple control:
    the gate is on while the phase (360 f t + ref) mod 360 lies in [alpha, alpha + width] degrees."""

    def __init__(self, anode, cathode, alpha, frequency, ref=0.0, width=120.0, Ron=1e-3, Roff=1e6):
        super().__init__(anode, cathode, Ron, Roff, 0.0)
        self.type = "THY"
        self.alpha, self.f, self.ref, self.width = alpha, frequency, ref, width

    def gate(self, t):
        theta = math.fmod(math.fmod(360.0*self.f*t + self.ref, 360.0) + 360.0, 360.0)
        return math.fmod(theta - self.alpha + 720.0, 360.0) <= self.width


# ---------------------------------------------------------------- machines

class SynchronousMachine(Element):
    """Synchronous machine, classical model: EMF E behind subtransient reactance xd'' (and ra), rotor
    angle and speed from the swing equation 2H dw/dt = Pm - Pe - D (w - 1) (per unit).
    nodes = a, b, c (terminals), xa, xb, xc (internal EMF nodes)."""

    state_vars = ("rl.vprev", "rl.iprev", "w", "delta", "Pe")

    def __init__(self, terminals, internal, S, V, f, H, xd2, ra, E, delta0, Pm, D=0.0):
        super().__init__("SM", list(terminals) + list(internal))
        if not (S > 0 and V > 0 and H > 0 and xd2 > 0):
            raise ValueError("SM: need S, V, H, xd'' > 0")
        self.S, self.ws, self.Vpk = S, 2*math.pi*f, V*math.sqrt(2.0/3.0)
        self.H, self.E, self.Pm, self.D = H, E, Pm, D
        Zb = V*V/S
        self.rl = RLBlock(np.eye(3)*ra*Zb, np.eye(3)*xd2*Zb/self.ws)
        self.tpos, self.xpos = [0, 1, 2], [3, 4, 5]
        self.w = 1.0
        self.delta = delta0*math.pi/180.0
        self.Pe = 0.0

    def stamp(self):
        g = np.zeros((6, 6))
        self.rl.stamp(g, self.xpos, self.tpos)
        return g

    def history(self, t):
        h = np.zeros(6)
        self.rl.history(h, self.xpos, self.tpos)
        return h

    def known_terminals(self):
        return self.xpos

    def known_voltages(self, t):
        th = self.ws*t + self.delta
        return [self.E*self.Vpk*math.cos(th - 2*math.pi*k/3) for k in range(3)]

    def update(self, t, v):
        self.rl.update(self.rl.branch_voltage(v, self.xpos, self.tpos))   # current from the EMF to the terminal
        e = self.known_voltages(t)
        self.Pe = (e[0]*self.rl.iprev[0] + e[1]*self.rl.iprev[1] + e[2]*self.rl.iprev[2])/self.S
        self.w += config.step*(self.Pm - self.Pe - self.D*(self.w - 1.0))/(2.0*self.H)   # swing equation
        self.delta += config.step*self.ws*(self.w - 1.0)

    def current_names(self):
        n = self.name
        return [n + ".a", n + ".b", n + ".c", n + ".speed[pu]", n + ".delta[deg]", n + ".Pe[pu]"]

    def currents(self):
        return list(self.rl.iprev) + [self.w, self.delta*180.0/math.pi, self.Pe]


class InductionMachine(Element):
    """Induction machine (motor convention), voltage behind transient inductance in the stator frame:
        v = Rs i + L' di/dt + e',  e' = (Lm/Lr) d(psi_r)/dt,
        d(psi_r)/dt = -(Rr/Lr)(psi_r - Lm i) + j w_r psi_r,  J dw_m/dt = Te - T_load (T_load ~ w^2).
    Parameters in per unit on the machine rating. nodes = a, b, c, xa, xb, xc (internal)."""

    state_vars = ("rl.vprev", "rl.iprev", "wm", "psir", "isprev", "eprime", "Te")

    def __init__(self, terminals, internal, S, V, f, poles, rs, xls, rr, xlr, xm, H, Tload, speed0=0.0):
        super().__init__("IM", list(terminals) + list(internal))
        if not (S > 0 and V > 0 and poles >= 2 and poles % 2 == 0 and xm > 0 and H > 0):
            raise ValueError("IM: need S, V, H, xm > 0 and an even number of poles")
        ws, Zb = 2*math.pi*f, V*V/S
        self.P2 = poles/2.0
        self.wms = ws/self.P2
        Lls, Llr = xls*Zb/ws, xlr*Zb/ws
        self.Lm = xm*Zb/ws
        self.Lr = Llr + self.Lm
        self.Rr = rr*Zb
        Lp = Lls + self.Lm - self.Lm*self.Lm/self.Lr          # transient inductance L'
        self.rl = RLBlock(np.eye(3)*rs*Zb, np.eye(3)*Lp)
        self.J = 2.0*H*S/(self.wms*self.wms)
        self.Tb = S/self.wms
        self.Tload = Tload
        self.wm = float(speed0*self.wms)
        self.psir = 0j
        self.isprev = 0j
        self.eprime = 0j
        self.Te = 0.0
        self.tpos, self.xpos = [0, 1, 2], [3, 4, 5]

    def stamp(self):
        g = np.zeros((6, 6))
        self.rl.stamp(g, self.tpos, self.xpos)
        return g

    def history(self, t):
        h = np.zeros(6)
        self.rl.history(h, self.tpos, self.xpos)
        return h

    def known_terminals(self):
        return self.xpos

    def known_voltages(self, t):
        a2, a1 = cmath.rect(1.0, -2*math.pi/3), cmath.rect(1.0, 2*math.pi/3)
        return [self.eprime.real, (self.eprime*a2).real, (self.eprime*a1).real]

    def update(self, t, v):
        self.rl.update(self.rl.branch_voltage(v, self.tpos, self.xpos))   # current from the terminal into the machine
        a = cmath.rect(1.0, 2*math.pi/3)
        i = self.rl.iprev
        is_ = (2.0/3.0)*(i[0] + a*i[1] + a*a*i[2])
        h = config.step
        A = complex(-self.Rr/self.Lr, self.P2*self.wm)
        B = self.Rr*self.Lm/self.Lr
        psi = (self.psir*(1.0 + A*h/2.0) + B*h/2.0*(is_ + self.isprev))/(1.0 - A*h/2.0)
        self.Te = 1.5*self.P2*(self.Lm/self.Lr)*(psi.real*is_.imag - psi.imag*is_.real)
        Tl = self.Tload*self.Tb*(self.wm/self.wms)**2
        self.wm += h*(self.Te - Tl)/self.J
        A2 = complex(-self.Rr/self.Lr, self.P2*self.wm)
        self.eprime = (self.Lm/self.Lr)*(A2*psi + B*is_)
        self.psir = psi
        self.isprev = is_

    def current_names(self):
        n = self.name
        return [n + ".a", n + ".b", n + ".c", n + ".speed[pu]", n + ".torque[pu]"]

    def currents(self):
        return list(self.rl.iprev) + [self.wm/self.wms, self.Te/self.Tb]
