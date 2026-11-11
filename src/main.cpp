#include <chrono>
#include <exception>
#include <iostream>

#include "circuit.h"
#include "config.h"
#include "netlist.h"
#include "network.h"

// Re-solve a step at most this many times while non-linear elements change state
constexpr int MAX_ITERATIONS = 50;

// Usage:  emtsim                         -> circuit in circuit.cpp
//         emtsim circuits/grid100.net    -> circuit from a netlist file
//         emtsim circuits/grid100.net --chatter=off --interp=off
//                                        -> override the netlist's .options (chatter removal, interpolation)
int main(int argc, char** argv) {
    try {
        std::string path;
        int chatter = -1, interp = -1;   // -1: as the netlist says
        for (int k = 1; k < argc; ++k) {
            const std::string arg = argv[k];
            const std::size_t eq = arg.find('=');
            const std::string key = arg.substr(0, eq);
            const int on = eq == std::string::npos ? -1 : parseOnOff(arg.substr(eq + 1));
            if ((key == "--chatter" || key == "--interp") && on >= 0)
                (key == "--chatter" ? chatter : interp) = on;
            else if (!arg.empty() && arg[0] == '-')
                throw std::runtime_error("usage: emtsim [netlist.net] [--chatter=on|off] [--interp=on|off]");
            else
                path = arg;
        }
        std::cout << "\nStart running the simulation..... \n\n";
        double Time = config::srtTime;

        const auto wallStart = std::chrono::steady_clock::now();

        // Initialize the network with the circuit from a netlist file or circuit.cpp
        Circuit ckt = !path.empty() ? loadNetlist(path) : buildCircuit();
        if (chatter >= 0) config::chatter = chatter;
        if (interp >= 0) config::interp = interp;
        std::cout << "Chatter removal " << (config::chatter ? "on" : "off") << ", switching interpolation "
                  << (config::interp ? "on" : "off") << "\n\n";
        Network EMTDC(std::move(ckt.elements), ckt.userNodes);
        EMTDC.updateSwitches(Time);   // initial switch/fault states
        EMTDC.refactor();

        if (config::chatter || config::interp) {
            EMTDC.runWithEvents(MAX_ITERATIONS);
        } else {
            // plain trapezoidal loop; half a step of tolerance because Time accumulates round-off
            while (Time < config::stpTime - 0.5 * config::Dt) {
                Time += config::Dt;
                if (EMTDC.updateSwitches(Time)) EMTDC.refactor();
                EMTDC.solveStep(Time, Network::AllValves, MAX_ITERATIONS);
                EMTDC.recordv(Time);
            }
        }
        const int unconverged = EMTDC.unconverged;

        // Setup + time loop, excluding writing the results
        const double simMs =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - wallStart).count();

        EMTDC.writeResults("results.csv");
        EMTDC.writeCurrents("currents.csv");
        EMTDC.writeEvents("events.csv");
        std::cout << "\nDone: " << EMTDC.time().size() << " steps, " << EMTDC.numNodes() << " nodes";
        if (EMTDC.totalNodes() > EMTDC.numNodes())
            std::cout << " (+" << EMTDC.totalNodes() - EMTDC.numNodes() << " internal)";
        std::cout << ".\n"
                  << "Sparse LDL: G_UU has " << EMTDC.nnzG_UU() << " non-zeros, L has " << EMTDC.nnzFactor()
                  << " below the diagonal, " << EMTDC.factorisations() << " factorisations";
        if (EMTDC.resolves()) std::cout << " (" << EMTDC.resolves() << " re-solves for non-linear elements)";
        std::cout << "\n";
        if (unconverged)
            std::cout << "Warning: non-linear elements did not settle within " << MAX_ITERATIONS << " iterations in "
                      << unconverged << " steps\n";
        std::cout << "Simulation time: " << simMs << " ms\n"
                  << "Results written to results.csv, currents.csv and events.csv. Plot them with:  python plot.py\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
}
