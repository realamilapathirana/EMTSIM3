#include "switching.h"

#include <algorithm>
#include <cmath>
#include <deque>

#include "config.h"
#include "elements.h"
#include "network.h"

void interpolateState(const Vector& x, const Vector& y, double al, Vector& out) {
    out.resize(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = x[i] + al * (y[i] - x[i]);
}

std::vector<Step> stepsAfterSwitching() {
    if (!config::chatter) return {};
    return {{0.5, true}, {0.5, true}};
}

// After every step, the earliest switching inside it is found (Network::detect, from the elements'
// eventAlpha), the state is interpolated back to that instant and the element switches there. The
// steps from stepsAfterSwitching() follow; the regular trapezoidal step after them is interpolated
// back onto the output time grid. Grid points passed during the special steps are recorded from the
// last two of them. The same algorithm as run_with_events() in switching.py.
void runWithEvents(Network& net, int maxIterations) {
    const double dt = config::Dt, t0 = config::srtTime, eps = 1e-6 * dt;
    const long nsteps = std::lround((config::stpTime - t0) / dt);
    double a = t0, tP = t0;      // time of the present solution, and of P
    long kRec = 0;
    Vector P = net.save(), S, tmp;   // last accepted point, the new one, scratch
    std::deque<Step> pending;        // special steps still to take after a switching
    int stall = 0;
    std::vector<Element*> hits;

    auto lerp = [&](const Vector& x, const Vector& y, double al) -> Vector& { interpolateState(x, y, al, tmp); return tmp; };
    auto recordTo = [&](double tS) {   // the grid points up to tS, from P and S
        while (kRec < nsteps) {
            const double tk = t0 + static_cast<double>(kRec + 1) * dt;
            if (tk > tS + eps) break;
            const double al = tS - tP <= 1e-12 * dt ? 1.0 : (tk - tP) / (tS - tP);
            const bool between = std::fabs(al - 1.0) > 1e-12;
            if (between) net.load(lerp(P, S, al));
            net.recordv(tk);
            if (between) net.load(S);
            ++kRec;
        }
    };
    auto accept = [&](double t) { recordTo(t); std::swap(P, S); tP = a = t; };

    while (kRec < nsteps) {
        const bool special = !pending.empty();
        const Step step = special ? pending.front() : Step{1.0, false};
        const double h = step.fraction * dt, b = a + h;
        config::step = h;
        config::euler = step.euler;
        for (Element* e : net.events()) e->beginStep(a, b);
        const bool fallback = stall > 20;   // valves switching back and forth at one instant: let check() settle them
        net.solveStep(b, fallback ? Network::AllValves : special ? Network::TurnOnOnly : Network::NoValves, maxIterations);
        const double best = fallback ? -1.0 : net.detect(a, h, config::interp, hits);
        if (fallback) stall = 0;
        S = net.save();
        if (best < 0.0) {
            if (special) {
                pending.pop_front();
                if (!pending.empty()) { std::swap(P, S); tP = a = b; }   // recorded after the last special step
                else accept(b);
                continue;
            }
            const double g = t0 + std::floor((b - t0) / dt + 1e-6) * dt;   // last grid point <= b
            if (a + eps < g && g < b - eps) {                                // off the grid: back onto it
                S.swap(lerp(P, S, (g - a) / dt));
                net.load(S);
                accept(g);
            } else {
                accept(std::fabs(g - b) <= eps ? g : b);
            }
            continue;
        }
        const double te = a + best * h;   // a switching inside the step
        if (best < 1.0 - 1e-12) {
            S.swap(lerp(P, S, best));
            net.load(S);
        }
        stall = best <= 1e-9 ? stall + 1 : 0;
        accept(te);
        bool changed = false;
        for (Element* e : hits) changed = e->applyEvent(te) || changed;
        if (changed) {
            net.refactor();
            const std::vector<Step> after = stepsAfterSwitching();
            pending.assign(after.begin(), after.end());
        }
    }
    config::step = dt;
    config::euler = false;
}
