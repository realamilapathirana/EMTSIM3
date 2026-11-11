"""Switching interpolation and chatter removal: the time loop used when either is on, and the two
pieces of the method that are meant to be changed (interpolate, steps_after_switching).
See docs/switching.md. The same algorithm as src/switching.cpp."""
import math

import config


def interpolate(x, y, al):
    """Interpolation: the network state at fraction al between two solution points x and y."""
    return x + al*(y - x)


def steps_after_switching():
    """Chatter removal: the steps to take right after a switching, each (fraction of Dt, backward
    Euler?). Critical damping adjustment is two backward-Euler half steps."""
    return [(0.5, True), (0.5, True)] if config.chatter else []


def run_with_events(net, max_iterations=50):
    """After every step, the earliest switching inside it is found (Network.detect, from the elements'
    event_alpha), the state is interpolated back to that instant and the element switches there. The
    steps from steps_after_switching() follow; the regular trapezoidal step after them is interpolated
    back onto the output time grid. Grid points passed during the special steps are recorded from
    the last two of them."""
    dt, t0 = config.Dt, config.srtTime
    nsteps = int(round((config.stpTime - t0)/dt))
    eps = 1e-6*dt
    a = tP = t0            # time of the present solution, and of P
    k_rec = 0
    P = net.save()         # last accepted point; S is the new one
    S = None
    pending = []           # special steps still to take after a switching
    stall = 0

    def record_to(tS):     # the grid points up to tS, from P and S
        nonlocal k_rec
        while k_rec < nsteps:
            tk = t0 + (k_rec + 1)*dt
            if tk > tS + eps:
                break
            al = 1.0 if tS - tP <= 1e-12*dt else (tk - tP)/(tS - tP)
            between = abs(al - 1.0) > 1e-12
            if between:
                net.load(interpolate(P, S, al))
            net.recordv(tk)
            if between:
                net.load(S)
            k_rec += 1

    def accept(t):
        nonlocal P, tP, a
        record_to(t)
        P, tP, a = S, t, t

    while k_rec < nsteps:
        special = bool(pending)
        fraction, euler = pending[0] if special else (1.0, False)
        h = fraction*dt
        b = a + h
        config.step, config.euler = h, euler
        for obj in net.events:
            obj.begin_step(a, b)
        fallback = stall > 20     # valves switching back and forth at one instant: let check() settle them
        net.solve_step(b, valves=True if fallback else "on" if special else False, max_iterations=max_iterations)
        best, hits = (None, []) if fallback else net.detect(a, h, config.interp)
        if fallback:
            stall = 0
        S = net.save()
        if best is None:
            if special:
                pending.pop(0)
                if pending:        # recorded after the last special step
                    P, tP, a = S, b, b
                else:
                    accept(b)
                continue
            g = t0 + math.floor((b - t0)/dt + 1e-6)*dt   # last grid point <= b
            if a + eps < g < b - eps:                     # off the grid: back onto it
                S = interpolate(P, S, (g - a)/dt)
                net.load(S)
                accept(g)
            else:
                accept(g if abs(g - b) <= eps else b)
            continue
        te = a + best*h           # a switching inside the step
        if best < 1.0 - 1e-12:
            S = interpolate(P, S, best)
            net.load(S)
        stall = stall + 1 if best <= 1e-9 else 0
        accept(te)
        changed = False
        for obj in hits:
            changed = obj.apply_event(te) or changed
        if changed:
            net.refactor()
            pending = list(steps_after_switching())
    config.step, config.euler = dt, False
