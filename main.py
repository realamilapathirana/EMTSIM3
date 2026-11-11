import sys
import emtcls
import config
import elements as el
import dataloader as dl
from switching import run_with_events


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
        run_with_events(EMTDC, MAX_ITERATIONS)
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
