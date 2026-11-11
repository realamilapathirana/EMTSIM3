"""Read an EMTSim netlist file into a list of elements (same format as src/netlist.cpp).

One element per line; node 0 = ground; '#' starts a comment, and the first word of
the comment (if any) names the element in the recorded currents. Optional
parameters can be given as key=value.

    R  a b  ohm                     L  a b  henry                 C  a b  farad
    V  a b  peak angle freq [R] [L]     R > 0: Norton; L: ideal source behind series R + L
    SW a b  closed|open t1 t2 ...  [czero] [chop=A] [ron=] [roff=]
    F  a b  t_apply [t_clear]  [ron=] [roff=]
    RL3 a1 b1 c1 a2 b2 c2  R1 L1 [R0 L0]            coupled three-phase series R-L (sequence data)
    PI  a b  R L C                                   single-phase nominal-pi line (C total)
    PI3 a1 b1 c1 a2 b2 c2  R1 L1 C1 [R0 L0 C0]       three-phase nominal-pi line
    V3  a b c  peak angle freq [R1 L1 [R0 L0]]       three-phase source (behind impedance if given)
    TL  k m  Zc tau [R]                              Bergeron line
    TL3 ka kb kc ma mb mc  Z1 tau1 Z0 tau0 [R1 R0]   three-phase Bergeron line (transposed)
    XF  p q s t  S V1 V2 f x r  [im=%] [knee=pu lsat=pu]           single-phase transformer
    XF3 a1 b1 c1 a2 b2 c2  conn S V1 V2 f x r  [im=%] [knee= lsat=]  conn: YNyn0 YNd1 YNd11 Dyn1 Dyn11 Yy0 Dd0 ...
    LSAT a b  i1 flux1 i2 flux2 ...                  saturable inductor
    MOV a b  Vref [Iref] [q]                         surge arrester
    D   a k  [ron=] [roff=] [vf=]                    diode
    THY a k  alpha freq [ref=deg] [width=deg] [ron=] [roff=]   thyristor
    SM  a b c  S V f H xd2 ra E delta0 Pm [D=]       synchronous machine (classical)
    IM  a b c  S V f poles rs xls rr xlr xm H Tload [speed0=]   induction machine
    .tran dt tstop                                   time step and stop time (before the elements)
    .options chatter=on|off interp=on|off            chatter removal / switching interpolation (default on)
"""
import math
import re


import config
import elements as el

TERMINALS = {"R": 2, "L": 2, "C": 2, "V": 2, "SW": 2, "F": 2, "PI": 2, "TL": 2, "LSAT": 2, "MOV": 2, "D": 2,
             "THY": 2, "XF": 4, "RL3": 6, "PI3": 6, "TL3": 6, "XF3": 6, "V3": 3, "SM": 3, "IM": 3}


class ElementList(list):
    """List of elements; user_nodes = highest node number in the file (higher = internal nodes)"""
    user_nodes = None


def _num(s, what):
    try:
        return float(s)
    except ValueError:
        raise ValueError('%s "%s" is not a number' % (what, s))


class _Args:
    def __init__(self, args, kw):
        self.args, self.kw = args, kw

    def req(self, k, what):
        if k >= len(self.args):
            raise ValueError("missing " + what)
        return _num(self.args[k], what)

    def opt(self, k, key, default):
        if key in self.kw:
            return self.kw[key]
        return _num(self.args[k], key) if k < len(self.args) else default

    def kwor(self, key, default):
        return self.kw.get(key, default)

    def has(self, k):
        return k < len(self.args)


class _Builder:
    def __init__(self, max_node):
        self.elements = ElementList()
        self.elements.user_nodes = max_node
        self.next = max_node + 1
        self.used = set()

    def node(self):
        n = self.next
        self.next += 1
        return n

    def add(self, name, element, type_=None):
        n, k = name, 2
        while n in self.used:
            n = "%s~%d" % (name, k)
            k += 1
        self.used.add(n)
        element.name = n
        if type_:
            element.type = type_
        self.elements.append(element)
        return element


def _build(kind, name, n, args, kw, flags, b):
    a = _Args(args, kw)
    if kind in ("R", "L", "C"):
        b.add(name, el.RLC(brnType=kind, strnode=n[0], stpnode=n[1], value=a.req(0, "value")))
    elif kind == "V":
        mag, ang, f = a.req(0, "magnitude"), a.req(1, "angle"), a.req(2, "frequency")
        R, L = a.opt(3, "r", 0.0), a.opt(4, "l", 0.0)
        if L > 0:   # ideal source behind series R + L
            if n[0] != 0 and n[1] != 0:
                raise ValueError("a source with inductance needs one terminal on ground")
            x = b.node()
            if n[0] == 0:
                b.add(name, el.Source(strnode=0, stpnode=x, magnitude=mag, angle=ang, frequency=f))
                b.add(name + ".Z", el.CoupledRL([x], [n[1]], [[R]], [[L]]))
            else:
                b.add(name, el.Source(strnode=x, stpnode=0, magnitude=mag, angle=ang, frequency=f))
                b.add(name + ".Z", el.CoupledRL([n[0]], [x], [[R]], [[L]]))
        else:
            b.add(name, el.Source(strnode=n[0], stpnode=n[1], magnitude=mag, angle=ang, frequency=f, resistance=R))
    elif kind == "SW":
        if not args or args[0] not in ("closed", "open"):
            raise ValueError("SW state must be 'closed' or 'open'")
        b.add(name, el.Switch(strnode=n[0], stpnode=n[1], Ron=a.kwor("ron", 0.01), Roff=a.kwor("roff", 1.0e6),
                              closed=(args[0] == "closed"), t_switch=[_num(x, "switching time") for x in args[1:]],
                              opens_at_current_zero=("czero" in flags or "chop" in kw), chop=a.kwor("chop", 0.0)))
    elif kind == "F":
        t1 = a.req(1, "t_clear") if a.has(1) else None
        b.add(name, el.Fault(strnode=n[0], stpnode=n[1], Ron=a.kwor("ron", 0.01), Roff=a.kwor("roff", 1.0e6),
                             t_apply=a.req(0, "t_apply"), t_clear=t1))
    elif kind == "RL3":
        R1, L1 = a.req(0, "R1"), a.req(1, "L1")
        b.add(name, el.CoupledRL.from_sequence(n[0:3], n[3:6], R1, L1, a.opt(2, "r0", R1), a.opt(3, "l0", L1)), "RL3")
    elif kind == "PI":
        R, L, C = a.req(0, "R"), a.req(1, "L"), a.req(2, "C")
        b.add(name, el.CoupledRL([n[0]], [n[1]], [[R]], [[L]]), "PI")
        if C > 0:
            b.add(name + ".Ck", el.RLC(brnType="C", strnode=n[0], stpnode=0, value=C/2.0))
            b.add(name + ".Cm", el.RLC(brnType="C", strnode=n[1], stpnode=0, value=C/2.0))
    elif kind == "PI3":
        R1, L1, C1 = a.req(0, "R1"), a.req(1, "L1"), a.req(2, "C1")
        R0, L0, C0 = a.opt(3, "r0", R1), a.opt(4, "l0", L1), a.opt(5, "c0", C1)
        b.add(name, el.CoupledRL.from_sequence(n[0:3], n[3:6], R1, L1, R0, L0), "PI3")
        if C1 > 0:
            b.add(name + ".Ck", el.CoupledC.from_sequence(n[0:3], C1/2.0, C0/2.0))
            b.add(name + ".Cm", el.CoupledC.from_sequence(n[3:6], C1/2.0, C0/2.0))
    elif kind == "V3":
        mag, ang, f = a.req(0, "magnitude"), a.req(1, "angle"), a.req(2, "frequency")
        shift, ph = (0.0, -120.0, 120.0), (".a", ".b", ".c")
        has_l = a.has(4) or "l1" in kw            # an inductance is given
        R1, L1 = a.opt(3, "r1", 0.0), a.opt(4, "l1", 0.0)
        if has_l:   # ideal EMFs behind a coupled (sequence) impedance
            x = []
            for p in range(3):
                x.append(b.node())
                b.add(name + ph[p], el.Source(strnode=0, stpnode=x[p], magnitude=mag, angle=ang + shift[p], frequency=f))
            b.add(name + ".Z", el.CoupledRL.from_sequence(x, n[0:3], R1, L1, a.opt(5, "r0", R1), a.opt(6, "l0", L1)))
        else:
            for p in range(3):
                b.add(name + ph[p], el.Source(strnode=0, stpnode=n[p], magnitude=mag, angle=ang + shift[p], frequency=f,
                                              resistance=R1))
    elif kind == "TL":
        b.add(name, el.BergeronLine([n[0]], [n[1]], [a.req(0, "Zc")], [a.req(1, "tau")], [a.opt(2, "r", 0.0)]))
    elif kind == "TL3":
        Z1, t1, Z0, t0 = a.req(0, "Z1"), a.req(1, "tau1"), a.req(2, "Z0"), a.req(3, "tau0")
        R1 = a.opt(4, "r1", 0.0)
        R0 = a.opt(5, "r0", R1)
        b.add(name, el.BergeronLine(n[0:3], n[3:6], [Z0, Z1, Z1], [t0, t1, t1], [R0, R1, R1]), "TL3")
    elif kind in ("XF", "XF3"):
        three = kind == "XF3"
        o = 1 if three else 0
        if three and not args:
            raise ValueError("XF3 needs a connection, e.g. YNd1")
        S, V1, V2, f = a.req(o + 0, "S"), a.req(o + 1, "V1"), a.req(o + 2, "V2"), a.req(o + 3, "f")
        x, r = a.req(o + 4, "x"), a.req(o + 5, "r")
        im, knee, lsat = a.kwor("im", 1.0), a.kwor("knee", 0.0), a.kwor("lsat", 2.0*x)
        w = 2.0*math.pi*f
        windings, delta_nodes = [], []
        Vw1, Vw2, Sw = V1, V2, S
        if not three:
            windings.append((n[0], n[1], n[2], n[3]))
        else:
            m = re.fullmatch(r"(YN|Y|D)(yn|y|d)(0|1|11)", args[0])
            if not m:
                raise ValueError('XF3 connection "%s": use e.g. YNyn0, YNd1, YNd11, Dyn1, Dyn11, Dd0' % args[0])
            P, Q, clock = m.group(1), m.group(2), int(m.group(3))
            pd, sd = P == "D", Q == "d"
            if (pd == sd and clock != 0) or (pd != sd and clock == 0):
                raise ValueError('XF3 connection "%s": star-star and delta-delta need clock 0, star-delta 1 or 11' % args[0])
            Vw1 = V1 if pd else V1/math.sqrt(3.0)
            Vw2 = V2 if sd else V2/math.sqrt(3.0)
            Sw = S/3.0
            Np = b.node() if P == "Y" else 0
            Ns = b.node() if Q == "y" else 0
            for p in range(3):
                nxt, prv = (p + 1) % 3, (p + 2) % 3
                # delta winding orientation sets the phase shift (clock 1: secondary lags 30 deg, 11: leads)
                pt = n[nxt if (sd or clock == 11) else prv] if pd else Np
                st = n[3 + (nxt if (pd or clock == 1) else prv)] if sd else Ns
                windings.append((n[p], pt, n[3 + p], st))
            # a delta winding has no ground reference of its own: add a very high resistance to ground
            if pd:
                delta_nodes += n[0:3]
            if sd:
                delta_nodes += n[3:6]
        Zb, ratio = Vw1*Vw1/Sw, Vw1/Vw2
        for p, wd in enumerate(windings):
            nm = name + (".a", ".b", ".c")[p] if three else name
            b.add(nm, el.TransformerCore(wd[0], wd[1], wd[2], wd[3], ratio, r*Zb, x*Zb/w), kind)
            if im > 0:
                Lm = Vw1*Vw1/(w*Sw*im/100.0)
                if knee > 0:   # two-slope saturation curve: Lm up to the knee, then lsat
                    lamR = math.sqrt(2.0)*Vw1/w
                    iK = knee*lamR/Lm
                    Ls = lsat*Zb/w
                    b.add(nm + ".m", el.SaturableInductor(wd[0], wd[1], [iK, iK + lamR/Ls], [knee*lamR, (knee + 1.0)*lamR]))
                else:
                    b.add(nm + ".m", el.RLC(brnType="L", strnode=wd[0], stpnode=wd[1], value=Lm))
        for d in delta_nodes:
            b.add(name + ".ref", el.RLC(brnType="R", strnode=d, stpnode=0, value=1.0e8))
    elif kind == "LSAT":
        if len(args) < 2 or len(args) % 2:
            raise ValueError("LSAT needs pairs of current and flux")
        b.add(name, el.SaturableInductor(n[0], n[1], [_num(v, "current") for v in args[0::2]],
                                         [_num(v, "flux") for v in args[1::2]]))
    elif kind == "MOV":
        b.add(name, el.SurgeArrester(n[0], n[1], a.req(0, "Vref"), a.opt(1, "iref", 1000.0), a.opt(2, "q", 26.0)))
    elif kind == "D":
        b.add(name, el.Diode(n[0], n[1], a.kwor("ron", 1e-3), a.kwor("roff", 1e6), a.kwor("vf", 0.0)))
    elif kind == "THY":
        b.add(name, el.Thyristor(n[0], n[1], a.req(0, "alpha"), a.req(1, "frequency"), a.kwor("ref", 0.0),
                                 a.kwor("width", 120.0), a.kwor("ron", 1e-3), a.kwor("roff", 1e6)))
    elif kind == "SM":
        x = [b.node(), b.node(), b.node()]
        b.add(name, el.SynchronousMachine(n[0:3], x, a.req(0, "S"), a.req(1, "V"), a.req(2, "f"), a.req(3, "H"),
                                          a.req(4, "xd2"), a.req(5, "ra"), a.req(6, "E"), a.req(7, "delta0"),
                                          a.req(8, "Pm"), a.opt(9, "d", 0.0)))
    elif kind == "IM":
        x = [b.node(), b.node(), b.node()]
        b.add(name, el.InductionMachine(n[0:3], x, a.req(0, "S"), a.req(1, "V"), a.req(2, "f"), int(a.req(3, "poles")),
                                        a.req(4, "rs"), a.req(5, "xls"), a.req(6, "rr"), a.req(7, "xlr"), a.req(8, "xm"),
                                        a.req(9, "H"), a.req(10, "Tload"), a.opt(11, "speed0", 0.0)))
    else:
        raise ValueError("unknown element type '%s'" % kind)


def parse_on_off(value):
    """on/1/true/yes -> True, off/0/false/no -> False, otherwise None"""
    v = value.lower()
    return True if v in ("on", "1", "true", "yes") else False if v in ("off", "0", "false", "no") else None


def parse(text, where="netlist"):
    statements, max_node = [], 0
    # time step and stop time: defaults unless the netlist has a ".tran dt tstop" line
    config.Dt = config.step = config.defaultDt
    config.euler = False
    config.stpTime = config.srtTime + config.defaultStop
    config.chatter, config.interp = config.defaultChatter, config.defaultInterp
    for lineno, raw in enumerate(text.splitlines(), 1):
        loc = "%s:%d" % (where, lineno)
        body, _, comment = raw.partition("#")
        tok = body.split()
        if not tok:
            continue
        if tok[0].upper() == ".TRAN":
            if len(tok) != 3:
                raise ValueError("%s: use .tran dt tstop" % loc)
            dt, stop = _num(tok[1], "dt"), _num(tok[2], "tstop")
            if not (dt > 0 and stop > dt):
                raise ValueError("%s: .tran needs 0 < dt < tstop" % loc)
            if statements:
                raise ValueError("%s: put .tran before the elements" % loc)
            config.Dt = config.step = dt
            config.stpTime = config.srtTime + stop
            continue
        if tok[0].upper() == ".OPTIONS":
            for kv in tok[1:]:
                key, _, val = kv.partition("=")
                on = parse_on_off(val)
                if key.lower() not in ("chatter", "interp") or on is None:
                    raise ValueError("%s: use .options chatter=on|off interp=on|off" % loc)
                setattr(config, key.lower(), on)
            continue
        cw = comment.split()
        name = cw[0].rstrip(",;:") if cw else ""
        kind = tok[0].upper()
        nt = TERMINALS.get(kind)
        try:
            if nt is None:
                raise ValueError("unknown element type '%s'" % tok[0])
            if len(tok) < 1 + nt:
                raise ValueError("expected %d node numbers" % nt)
            nodes = []
            for s in tok[1:1 + nt]:
                v = _num(s, "node")
                if v < 0 or v != int(v):
                    raise ValueError("node numbers must be whole numbers >= 0")
                nodes.append(int(v))
                max_node = max(max_node, int(v))
            args, kw, flags = [], {}, set()
            for s in tok[1 + nt:]:
                if "=" in s:
                    key, _, val = s.partition("=")
                    kw[key.lower()] = _num(val, key)
                elif s == "czero":
                    flags.add(s)
                else:
                    args.append(s)
        except ValueError as e:
            raise ValueError("%s: %s" % (loc, e))
        statements.append((kind, name or "%s@%d" % (kind, lineno), nodes, args, kw, flags, loc))

    b = _Builder(max_node)
    for kind, name, nodes, args, kw, flags, loc in statements:
        try:
            _build(kind, name, nodes, args, kw, flags, b)
        except ValueError as e:
            raise ValueError("%s: %s" % (loc, e))
    return b.elements


def load(path):
    with open(path) as f:
        return parse(f.read(), path)
