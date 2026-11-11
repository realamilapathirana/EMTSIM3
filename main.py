import math
import sys
import emtcls
import config
import elements as el
import dataloader as dl


# Re-solve a step at most this many times while non-linear elements change state
MAX_ITERATIONS = 50


def run(netlist_path=None, chatter=None, interp=None):
    """chatter / interp: True or False override the netlist's .options line (None: keep it)"""
    print("\nStart running the simulation..... \n")
    Time=config.srtTime
    #Use a netlist file if given, otherwise the circuit in dataloader.py
    if netlist_path:
        import netlist
        dl.elemet_list = netlist.load(netlist_path)
    if chatter is not None:
        config.chatter = chatter
    if interp is not None:
        config.interp = interp
    print("Chatter removal %s, switching interpolation %s\n" % ("on" if config.chatter else "off", "on" if config.interp else "off"))
    EMTDC=emtcls.Network()
    EMTDC.updateSwitches(Time)      # initial switch/fault states
    EMTDC.refactor()

    if config.chatter or config.interp:
        run_with_events(EMTDC)
    else:
        # plain trapezoidal loop; half a step of tolerance because Time accumulates round-off
        while Time < config.stpTime - 0.5*config.Dt:
            Time += config.Dt
            if EMTDC.updateSwitches(Time):
                EMTDC.refactor()
            EMTDC.solve_step(Time, valves=True, max_iterations=MAX_ITERATIONS)
            EMTDC.recordv(Time)
    finish(EMTDC)


def finish(EMTDC):
    if EMTDC.unconverged:
        print("Warning: non-linear elements did not settle within %d iterations in %d steps"
              % (MAX_ITERATIONS, EMTDC.unconverged))
    EMTDC.plotv()


def run_with_events(EMTDC):
    """Time loop with switching interpolation and chatter removal (Watson & Arrillaga 9.4.1); the same
    algorithm as Network::runWithEvents in src/network.cpp. After a step, the earliest switching inside
    it is found, the state is interpolated back to that instant and the element switches there. After a
    switching the next two steps are backward-Euler half steps (same G, no voltage term in the histories),
    which stops the trapezoidal chatter; the step after them is interpolated back onto the regular grid."""
    dt, t0 = config.Dt, config.srtTime
    nsteps = int(round((config.stpTime - t0)/dt))
    eps = 1e-6*dt
    a = tP = t0            # time of the present solution, and of P
    k_rec = 0
    P = EMTDC.save()       # last accepted point; S is the new one
    S = None
    be_left = stall = 0

    def record_to(tS):     # the grid points up to tS, from P and S
        nonlocal k_rec
        while k_rec < nsteps:
            tk = t0 + (k_rec + 1)*dt
            if tk > tS + eps:
                break
            al = 1.0 if tS - tP <= 1e-12*dt else (tk - tP)/(tS - tP)
            between = abs(al - 1.0) > 1e-12
            if between:
                EMTDC.load(P + al*(S - P))
            EMTDC.recordv(tk)
            if between:
                EMTDC.load(S)
            k_rec += 1

    def accept(t):
        nonlocal P, tP, a
        record_to(t)
        P, tP, a = S, t, t

    while k_rec < nsteps:
        be = be_left > 0
        h = 0.5*dt if be else dt
        b = a + h
        config.step, config.euler = h, be
        for obj in EMTDC.events:
            obj.begin_step(a, b)
        fallback = stall > 20     # valves switching back and forth at one instant: let check() settle them
        EMTDC.solve_step(b, valves=True if fallback else "on" if be else False, max_iterations=MAX_ITERATIONS)
        best, hits = (None, []) if fallback else EMTDC.detect(a, h, config.interp)
        if fallback:
            stall = 0
        S = EMTDC.save()
        if best is None:
            if be:
                be_left -= 1
                if be_left:        # first half step: its grid points are recorded after the second
                    P, tP, a = S, b, b
                else:
                    accept(b)
                continue
            g = t0 + math.floor((b - t0)/dt + 1e-6)*dt   # last grid point <= b
            if a + eps < g < b - eps:                     # off the grid: back onto it
                S = P + ((g - a)/dt)*(S - P)
                EMTDC.load(S)
                accept(g)
            else:
                accept(g if abs(g - b) <= eps else b)
            continue
        te = a + best*h           # a switching inside the step
        if best < 1.0 - 1e-12:
            S = P + best*(S - P)
            EMTDC.load(S)
        stall = stall + 1 if best <= 1e-9 else 0
        accept(te)
        changed = False
        for obj in hits:
            changed = obj.apply_event(te) or changed
        if changed:
            EMTDC.refactor()
            if config.chatter:
                be_left = 2
    config.step, config.euler = dt, False


if __name__ == "__main__":
    # python main.py                       -> circuit in dataloader.py
    # python main.py circuits/grid100.net  -> circuit from a netlist file
    # python main.py circuits/grid100.net --chatter=off --interp=off
    #                                      -> override the netlist's .options (chatter removal, interpolation)
    from netlist import parse_on_off
    path, opts = None, {}
    for arg in sys.argv[1:]:
        key, _, val = arg.partition("=")
        on = parse_on_off(val)
        if key in ("--chatter", "--interp") and on is not None:
            opts[key[2:]] = on
        elif arg.startswith("-"):
            raise SystemExit("usage: python main.py [netlist.net] [--chatter=on|off] [--interp=on|off]")
        else:
            path = arg
    run(path, **opts)
