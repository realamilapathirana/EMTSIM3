// Checks the solver against analytical solutions of first-order circuits.
// Run with:  ctest --test-dir build   (or run build/test_emtsim directly)

#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "circuit.h"
#include "config.h"
#include "linalg.h"
#include "netlist.h"
#include "network.h"

#include <algorithm>
#include <fstream>

namespace {

using cd = std::complex<double>;
constexpr double PI = 3.14159265358979323846;
constexpr double Vm = 100.0, f = 50.0, a_deg = -90.0;  // E(0) = 0
const double w = 2 * PI * f;
const cd ph = std::polar(1.0, a_deg * PI / 180.0);
const cd j(0.0, 1.0);

double E(double t) { return Vm * std::cos(w * t + a_deg * PI / 180.0); }

// Voltage across L in a series RL circuit energised at t = 0
double rl_vL(double t, double R, double L) {
    const cd I = Vm * ph / (R + j * w * L);
    const double i = std::real(I * std::exp(j * w * t)) - std::real(I) * std::exp(-R * t / L);
    return E(t) - R * i;
}

// Voltage across C in a series RC circuit energised at t = 0
double rc_vC(double t, double R, double C) {
    const cd Zc = 1.0 / (j * w * C);
    const cd Vc = Vm * ph * Zc / (R + Zc);
    return std::real(Vc * std::exp(j * w * t)) - std::real(Vc) * std::exp(-t / (R * C));
}

struct Run {
    std::vector<double> t;
    std::vector<double> v;
};

// Same time loop as main.cpp, for a fixed number of steps; returns the network with its recordings
std::unique_ptr<Network> runNetwork(Circuit ckt, int steps = 1000) {
    std::ostringstream quiet;
    auto* old = std::cout.rdbuf(quiet.rdbuf());   // silence matrix printing
    auto net = std::make_unique<Network>(std::move(ckt.elements), ckt.userNodes);
    double t = config::srtTime;
    net->updateSwitches(t);
    net->refactor();
    for (int s = 0; s < steps; ++s) {
        t += config::Dt;
        if (net->updateSwitches(t)) net->refactor();
        net->solveStep(t, Network::AllValves, 50);
        net->recordv(t);
    }
    std::cout.rdbuf(old);
    return net;
}

// Recorded voltage of one node
Run simulate(Circuit ckt, int node, int steps = 1000) {
    const auto net = runNetwork(std::move(ckt), steps);
    Run r;
    r.t = net->time();
    for (const auto& row : net->trace()) r.v.push_back(node == 0 ? 0.0 : row[node - 1]);
    return r;
}

int failures = 0;

void check(bool ok, const std::string& what, double value) {
    std::printf("%s  %-58s %.2e\n", ok ? "PASS" : "FAIL", what.c_str(), value);
    if (!ok) ++failures;
}

double maxErr(const Run& r, const std::function<double(double)>& exact, double tmin = 0.0) {
    double e = 0.0;
    for (std::size_t i = 0; i < r.t.size(); ++i)
        if (r.t[i] > tmin) e = std::max(e, std::fabs(r.v[i] - exact(r.t[i])) / Vm);
    return e;
}

// ---------------------------------------------------------------- helpers for the element tests

// Run a netlist for `steps` steps
std::unique_ptr<Network> runText(const std::string& text, int steps = 1000) { return runNetwork(parseNetlist(text, "test"), steps); }

std::vector<double> nodeTrace(const Network& net, int node) {
    std::vector<double> v;
    for (const auto& row : net.trace()) v.push_back(node == 0 ? 0.0 : row[node - 1]);
    return v;
}
std::vector<double> currentTrace(const Network& net, const std::string& name) {
    const auto& names = net.currentNames();
    const auto it = std::find(names.begin(), names.end(), name);
    if (it == names.end()) throw std::runtime_error("no recorded current " + name);
    const std::size_t k = static_cast<std::size_t>(it - names.begin());
    std::vector<double> v;
    for (const auto& row : net.currentTrace()) v.push_back(row[k]);
    return v;
}
// Phasor (amplitude, phase in degrees) of x at frequency f, fitted over the last `cycles` cycles
cd phasorOf(const std::vector<double>& t, const std::vector<double>& x, double fr, int cycles = 1) {
    const double T = cycles / fr, tend = t.back();
    double cc = 0, ss = 0, cs = 0, xc = 0, xs = 0;
    for (std::size_t i = 0; i < t.size(); ++i) {
        if (t[i] <= tend - T + 1e-12) continue;
        const double c = std::cos(2 * PI * fr * t[i]), s = std::sin(2 * PI * fr * t[i]);
        cc += c * c; ss += s * s; cs += c * s; xc += x[i] * c; xs += x[i] * s;
    }
    const double det = cc * ss - cs * cs, A = (xc * ss - xs * cs) / det, B = (xs * cc - xc * cs) / det;
    return cd(A, -B);   // x = A cos(wt) + B sin(wt) = Re{(A - jB) e^{jwt}}
}
double peakAbs(const std::vector<double>& x, std::size_t from = 0, std::size_t to = SIZE_MAX) {
    double m = 0;
    for (std::size_t i = from; i < std::min(to, x.size()); ++i) m = std::max(m, std::fabs(x[i]));
    return m;
}
double relErr(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }
double angleDiff(double a, double b) { return std::remainder(a - b, 360.0); }

void newElementTests() {
    std::printf("\n-- new elements --\n");
    const double fr = 50.0;   // w = 2 pi fr is the file-level constant

    {   // branch currents: series R-L, KCL source current = resistor current, and the analytical current
        const auto net = runText("V 0 1 100 -90 50   # SRC\nR 1 2 2   # R1\nL 2 0 10e-3   # L1\n");
        const auto iR = currentTrace(*net, "R1"), iS = currentTrace(*net, "SRC"), iL = currentTrace(*net, "L1");
        double kcl = 0, err = 0;
        const cd I = 100.0 * std::polar(1.0, -PI / 2) / (2.0 + j * w * 10e-3);
        for (std::size_t s = 0; s < iR.size(); ++s) {
            kcl = std::max(kcl, std::max(std::fabs(iR[s] - iS[s]), std::fabs(iR[s] - iL[s])));
            const double t = net->time()[s];
            err = std::max(err, std::fabs(iR[s] - (std::real(I * std::exp(j * w * t)) - std::real(I) * std::exp(-2.0 * t / 10e-3))));
        }
        check(kcl < 1e-9, "branch currents: ideal source = R = L current (KCL, max diff A)", kcl);
        check(err / std::abs(I) < 1e-3, "branch current of R-L matches the analytical current (max err pu)", err / std::abs(I));
    }
    {   // coupled R-L with equal sequence values == separate R and L per phase
        const auto a = runText("V3 1 2 3 100 0 50\nRL3 1 2 3 4 5 6 1.0 0.01\nR 4 0 10\nR 5 0 10\nR 6 0 10\n");
        const auto b = runText("V3 1 2 3 100 0 50\nR 1 7 1.0\nL 7 4 0.01\nR 2 8 1.0\nL 8 5 0.01\nR 3 9 1.0\nL 9 6 0.01\nR 4 0 10\nR 5 0 10\nR 6 0 10\n");
        double d = 0;
        for (int nd = 4; nd <= 6; ++nd) {
            const auto x = nodeTrace(*a, nd), y = nodeTrace(*b, nd);
            for (std::size_t s = 0; s < x.size(); ++s) d = std::max(d, std::fabs(x[s] - y[s]));
        }
        check(d < 1e-9, "RL3 with R0 = R1, L0 = L1 equals separate R and L in each phase (max diff V)", d);
    }
    {   // single-line-to-ground fault through a coupled line: If = 3E / (2 Z1 + Z0); healthy phase rises
        const double R1 = 1, L1 = 0.01, R0 = 3, L0 = 0.03, E = 100;
        const auto net = runText("V3 1 2 3 100 0 50\nRL3 1 2 3 4 5 6 1 0.01 3 0.03\nF 4 0 0.0 ron=1e-6   # FLT\n", 1500);
        const cd Z1(R1, w * L1), Z0(R0, w * L0), If = 3.0 * E / (2.0 * Z1 + Z0);
        const cd Zm = (Z0 - Z1) / 3.0, Vb = E * std::polar(1.0, -2 * PI / 3) - Zm * If;
        const auto iF = currentTrace(*net, "FLT");
        const double amp = std::abs(phasorOf(net->time(), iF, fr));
        check(relErr(amp, std::abs(If)) < 2e-3, "RL3 single-line-to-ground fault current = 3E/(2Z1+Z0) (rel err)", relErr(amp, std::abs(If)));
        const double vb = std::abs(phasorOf(net->time(), nodeTrace(*net, 5), fr));
        check(relErr(vb, std::abs(Vb)) < 2e-3, "RL3: healthy phase b voltage includes the mutual coupling (rel err)", relErr(vb, std::abs(Vb)));
    }
    {   // three-phase source behind sequence impedance: same fault, same formula
        const auto net = runText("V3 1 2 3 100 0 50 1 0.01 3 0.03\nF 1 0 0.0 ron=1e-6   # FLT\n", 1500);
        const cd If = 3.0 * 100.0 / (2.0 * cd(1, w * 0.01) + cd(3, w * 0.03));
        const double amp = std::abs(phasorOf(net->time(), currentTrace(*net, "FLT"), fr));
        check(relErr(amp, std::abs(If)) < 2e-3, "V3 with source impedance: SLG fault current = 3E/(2Z1+Z0) (rel err)", relErr(amp, std::abs(If)));
    }
    {   // nominal-pi line, open at the far end: V_r = V_s / (1 + Z Y / 2)
        const double R = 2, L = 0.05, C = 40e-6;
        const auto net = runText("V 0 1 100 0 50\nPI 1 2 2 0.05 40e-6\nR 2 0 1e6\n", 2000);
        const cd Z(R, w * L), Y(0, w * C), Vr = 100.0 / (1.0 + Z * Y / 2.0);
        const double amp = std::abs(phasorOf(net->time(), nodeTrace(*net, 2), fr));
        check(relErr(amp, std::abs(Vr)) < 2e-3, "PI line open-circuit (Ferranti) voltage = Vs / (1 + ZY/2) (rel err)", relErr(amp, std::abs(Vr)));
        const auto n3 = runText("V3 1 2 3 100 0 50\nPI3 1 2 3 4 5 6 2 0.05 40e-6 6 0.15 25e-6\nR 4 0 1e6\nR 5 0 1e6\nR 6 0 1e6\n", 2000);
        const double a3 = std::abs(phasorOf(n3->time(), nodeTrace(*n3, 4), fr));
        check(relErr(a3, std::abs(Vr)) < 2e-3, "PI3 with balanced excitation behaves as the positive-sequence line (rel err)", relErr(a3, std::abs(Vr)));
    }
    {   // breaker opening at current zero
        const auto net = runText("V 0 1 100 0 50\nSW 1 2 closed 0.0305 czero   # BRK\nR 2 3 1\nL 3 0 0.01\n");
        const auto i = currentTrace(*net, "BRK");
        const auto& t = net->time();
        std::size_t open = 0;
        for (std::size_t s = 1; s < i.size(); ++s)
            if (t[s] > 0.0305 && std::fabs(i[s]) < 1e-3 * peakAbs(i) && std::fabs(i[s - 1]) < 0.05 * peakAbs(i)) { open = s; break; }
        const bool ok = open && t[open] > 0.0305 && t[open] < 0.0305 + 0.0101 && peakAbs(i, open) < 1e-3 * peakAbs(i);
        check(ok, "breaker with czero opens at the first current zero after the command (opening time s)", open ? t[open] : -1.0);
        const auto n2 = runText("V 0 1 100 0 50\nSW 1 2 closed 0.0305   # BRK\nR 2 3 1\nL 3 0 0.01\n");
        const auto i2 = currentTrace(*n2, "BRK");
        check(std::fabs(i2[303]) > 0.3 * peakAbs(i2) && std::fabs(i2[305]) < 1e-3 * peakAbs(i2),
              "without czero the same breaker interrupts a large current at the command time (A just before)", i2[303]);
    }
    {   // Bergeron line: matched source, open end -> v_m(t) = E(t - tau) exactly (lossless, integer delay)
        const double Zc = 400, tau = 1e-3;
        const auto net = runText("V 0 1 100 -90 50\nR 1 2 400\nTL 2 3 400 1e-3\nR 3 0 1e9\n");
        double err = 0;
        const auto v = nodeTrace(*net, 3);
        for (std::size_t s = 0; s < v.size(); ++s) {
            const double t = net->time()[s], ex = t >= tau ? 100.0 * std::cos(w * (t - tau) - PI / 2) : 0.0;
            err = std::max(err, std::fabs(v[s] - ex));
        }
        check(err < 0.05, "Bergeron line, open end: v_m(t) = E(t - tau) (max err V of 100 V)", err);
        (void)Zc;
        const auto n2 = runText("V 0 1 100 -90 50\nR 1 2 400\nTL 2 3 400 1.05e-3\nR 3 0 1e9\n");
        double e2 = 0;
        const auto v2 = nodeTrace(*n2, 3);
        for (std::size_t s = 0; s < v2.size(); ++s) {
            const double t = n2->time()[s], ex = t >= 1.05e-3 + 1e-4 ? 100.0 * std::cos(w * (t - 1.05e-3) - PI / 2) : 0.0;
            if (t >= 1.05e-3 + 1e-4) e2 = std::max(e2, std::fabs(v2[s] - ex));
        }
        check(e2 < 2.0, "Bergeron line with a non-integer delay (interpolated): max err V of 100 V", e2);
        // lossy line at DC: I = E / (Rs + R + Rload) after the reflections settle
        const auto n3 = runText("V 0 1 100 0 0\nR 1 2 50\nTL 2 3 400 1e-4 20   # LINE\nR 3 0 30   # LOAD\n", 3000);
        const double iend = currentTrace(*n3, "LOAD").back();
        check(relErr(iend, 100.0 / (50 + 20 + 30)) < 5e-3, "lossy Bergeron line at DC: I = E / (Rs + R + Rload) (rel err)", relErr(iend, 1.0));
        // TL3 with balanced excitation == single-phase line with Z1, tau1
        const auto a3 = runText("V3 1 2 3 100 -90 50\nR 1 4 300\nR 2 5 300\nR 3 6 300\nTL3 4 5 6 7 8 9 300 1e-3 600 1.5e-3\nR 7 0 1e9\nR 8 0 1e9\nR 9 0 1e9\n");
        const auto a1 = runText("V 0 1 100 -90 50\nR 1 2 300\nTL 2 3 300 1e-3\nR 3 0 1e9\n");
        double d = 0;
        const auto x = nodeTrace(*a3, 7), y = nodeTrace(*a1, 3);
        for (std::size_t s = 0; s < x.size(); ++s) d = std::max(d, std::fabs(x[s] - y[s]));
        check(d < 1e-6, "TL3 with balanced excitation equals the single-phase line (Z1, tau1) (max diff V)", d);
    }
    {   // transformer: ratio and short-circuit current
        const double S = 1e6, V1 = 11e3, V2 = 400, x = 0.06, r = 0.01;
        const std::string xf = " 1e6 11000 400 50 0.06 0.01";
        const auto oc = runText("V 0 1 15556.3 0 50\nXF 1 0 2 0" + xf + "\nR 2 0 1e6\n", 2000);
        const double v2 = std::abs(phasorOf(oc->time(), nodeTrace(*oc, 2), fr));
        check(relErr(v2, 15556.3 * V2 / V1) < 1e-3, "transformer open-circuit secondary voltage = V1 / n (rel err)", relErr(v2, 15556.3 * V2 / V1));
        const auto sc = runText("V 0 1 15556.3 0 50\nXF 1 0 2 0" + xf + " im=0   # T\nR 2 0 1e-9\n", 2000);
        const double Zb = V1 * V1 / S, isc = 15556.3 / std::abs(cd(r * Zb, x * Zb));
        const double ip = std::abs(phasorOf(sc->time(), currentTrace(*sc, "T.1"), fr));
        check(relErr(ip, isc) < 2e-3, "transformer short-circuit current = V / |Zleak| (rel err)", relErr(ip, isc));
        const double is = std::abs(phasorOf(sc->time(), currentTrace(*sc, "T.2"), fr));
        check(relErr(is, isc * V1 / V2) < 2e-3, "transformer secondary current = n x primary current (rel err)", relErr(is, isc * V1 / V2));
        (void)S;
    }
    {   // three-phase transformer vector groups: phase shift of secondary phase a against primary phase a
        struct Case { const char* conn; double shift; };
        const Case cases[] = {{"YNyn0", 0}, {"YNd1", -30}, {"YNd11", 30}, {"Dyn1", -30}, {"Dyn11", 30}, {"Dd0", 0}, {"Yy0", 0}};
        for (const auto& c : cases) {
            const auto net = runText(std::string("V3 1 2 3 8981.6 0 50\nXF3 1 2 3 4 5 6 ") + c.conn +
                                     " 1e6 11000 400 50 0.06 0.01\nR 4 0 1e4\nR 5 0 1e4\nR 6 0 1e4\n", 2000);
            const cd p = phasorOf(net->time(), nodeTrace(*net, 1), fr), s = phasorOf(net->time(), nodeTrace(*net, 4), fr);
            const double sh = angleDiff(std::arg(s) * 180 / PI, std::arg(p) * 180 / PI);
            const double ratio = std::abs(s) / std::abs(p);
            const bool ok = std::fabs(angleDiff(sh, c.shift)) < 0.5 && relErr(ratio, 400.0 / 11000.0) < 5e-3;
            check(ok, std::string("XF3 ") + c.conn + ": secondary phase shift (deg) and voltage ratio", sh);
        }
    }
    {   // saturation: energising at voltage zero drives the flux to ~2x rated -> inrush much larger than without saturation
        const std::string base = "V 0 1 15556.3 -90 50\nXF 1 0 2 0 1e6 11000 400 50 0.06 0.01 im=1";
        const auto lin = runText(base + "   # T\nR 2 0 1e6\n", 600);
        const auto sat = runText(base + " knee=1.15 lsat=0.2   # T\nR 2 0 1e6\n", 600);
        const double Imag = 0.01 * 1e6 / 11000 * std::sqrt(2.0);   // rated magnetising current, peak
        const double pl = peakAbs(currentTrace(*lin, "T.m")), ps = peakAbs(currentTrace(*sat, "T.m"));
        check(pl < 2.2 * Imag && ps > 20 * Imag, "transformer inrush: linear core ~2x rated magnetising current, saturating core > 20x (x rated)", ps / Imag);
    }
    {   // saturable inductor: current follows the piecewise-linear flux-current curve
        const auto net = runText("V 0 1 400 0 50\nR 1 2 0.01\nLSAT 2 0 1 1.0 10 1.5   # LS\n", 1000);
        const auto i = currentTrace(*net, "LS"), v = nodeTrace(*net, 2);
        double lam = 0, vprev = 0, err = 0;
        for (std::size_t s = 0; s < i.size(); ++s) {
            lam += 0.5 * config::Dt * (vprev + v[s]);
            vprev = v[s];
            const double a = std::fabs(lam), sgn = lam >= 0 ? 1 : -1;
            const double fcurve = a <= 1.0 ? a : a <= 1.5 ? 1.0 + (a - 1.0) / 0.5 * 9.0 : 10.0 + (a - 1.5) / (0.5 / 9.0);
            err = std::max(err, std::fabs(i[s] - sgn * fcurve));
        }
        check(err < 1e-6 * std::max(1.0, peakAbs(i)), "saturable inductor follows its flux-current curve (max err A)", err);
    }
    {   // surge arrester: clamps the voltage and follows i = Iref (v/Vref)^q (piecewise-linear)
        const auto net = runText("V 0 1 2000 0 50\nR 1 2 1\nMOV 2 0 1000 1000 26   # A\n");
        const auto v = nodeTrace(*net, 2), i = currentTrace(*net, "A");
        double worst = 0;
        for (std::size_t s = 0; s < v.size(); ++s)
            if (std::fabs(i[s]) > 10.0) worst = std::max(worst, relErr(std::fabs(i[s]), 1000 * std::pow(std::fabs(v[s]) / 1000, 26.0)));
        check(peakAbs(v) < 1080 && worst < 0.05, "surge arrester: voltage clamped near Vref, current on the characteristic (max rel err)", worst);
    }
    {   // diode: half-wave rectifier; full-wave bridge
        const auto net = runText("V 0 1 100 0 50\nD 1 2\nR 2 0 10\n");
        double err = 0;
        const auto v = nodeTrace(*net, 2);
        for (std::size_t s = 0; s < v.size(); ++s)
            err = std::max(err, std::fabs(v[s] - std::max(0.0, 100.0 * std::cos(w * net->time()[s]))));
        check(err < 0.05, "diode half-wave rectifier: v_load = max(E, 0) (max err V of 100 V)", err);
        const auto br = runText("V 0 1 100 0 50\nD 1 2\nD 3 1\nD 0 2\nD 3 0\nR 2 3 10\n");
        double e2 = 0;
        const auto vp = nodeTrace(*br, 2), vn = nodeTrace(*br, 3);
        for (std::size_t s = 0; s < vp.size(); ++s)
            e2 = std::max(e2, std::fabs((vp[s] - vn[s]) - std::fabs(100.0 * std::cos(w * br->time()[s]))));
        check(e2 < 0.1, "diode bridge: v_dc = |E| (max err V of 100 V)", e2);
    }
    {   // thyristor: controlled half-wave rectifier, average load voltage = Vm (1 + cos alpha) / (2 pi)
        const auto net = runText("V 0 1 100 -90 50\nTHY 1 2 60 50 ref=0\nR 2 0 10\n", 1000);
        const auto v = nodeTrace(*net, 2);
        double sum = 0;
        int n = 0;
        for (std::size_t s = 0; s < v.size(); ++s)
            if (net->time()[s] > 0.02 + 1e-9) { sum += v[s]; ++n; }
        const double avg = sum / n, ex = 100.0 * (1 + std::cos(PI / 3)) / (2 * PI);
        check(relErr(avg, ex) < 0.01, "thyristor rectifier (alpha = 60 deg): average voltage = Vm (1 + cos a) / 2 pi (rel err)", relErr(avg, ex));
    }
    {   // synchronous machine on an infinite bus: equilibrium at P = E V sin(delta) / X, swing frequency
        const double S = 100e6, V = 13.8e3, xd2 = 0.2, xline = 0.3, E = 1.1, Pm = 0.8, H = 1.0;
        const double Zb = V * V / S, Xt = xd2 + xline, d0 = std::asin(Pm * Xt / E) * 180 / PI;
        const double Lline = xline * Zb / w, Vpk = V * std::sqrt(2.0 / 3.0);
        char buf[400];
        std::snprintf(buf, sizeof buf,
                      "V3 4 5 6 %.6f 0 50\nRL3 1 2 3 4 5 6 0 %.9g\nSM 1 2 3 1e8 13800 50 %g %g 0 %g %.6f %g d=0   # G\n",
                      Vpk, Lline, H, xd2, E, d0 + 5.0, Pm);
        const auto net = runText(buf, 12000);
        const auto d = currentTrace(*net, "G.delta[deg]");
        // swing frequency from upward crossings of the mean by delta averaged over one 50 Hz cycle
        // (averaging removes the 50 Hz ripple left by energising the network from zero)
        const std::size_t win = 200;
        std::vector<double> ds;
        for (std::size_t s = win; s < d.size(); ++s) {
            double a = 0;
            for (std::size_t k = s - win; k < s; ++k) a += d[k];
            ds.push_back(a / win);
        }
        double mean = 0;
        for (double x : ds) mean += x;
        mean /= double(ds.size());
        std::vector<double> ups;
        for (std::size_t s = 1; s < ds.size(); ++s)
            if (ds[s - 1] < mean && ds[s] >= mean) ups.push_back(net->time()[s + win]);
        const double Pmax = E / Xt, fn = std::sqrt(w * Pmax * std::cos(d0 * PI / 180) / (2 * H)) / (2 * PI);
        const double fsim = ups.size() >= 2 ? (ups.size() - 1) / (ups.back() - ups.front()) : 0.0;
        check(relErr(fsim, fn) < 0.05, "synchronous machine: swing frequency = sqrt(ws Pmax cos d0 / 2H) / 2 pi (rel err)", relErr(fsim, fn));
        double dmean = 0;
        for (std::size_t s = d.size() / 2; s < d.size(); ++s) dmean += d[s];
        dmean /= double(d.size() - d.size() / 2);
        check(std::fabs(dmean - d0) < 1.0, "synchronous machine: rotor angle oscillates about the equilibrium angle (deg off)", dmean - d0);
    }
    {   // induction motor: direct-on-line start, final speed matches the equivalent circuit at the load torque
        const double S = 1e5, V = 400, rs = 0.01, xls = 0.1, rr = 0.012, xlr = 0.1, xm = 3.0, H = 0.3, Tl = 0.5;
        char buf[300];
        std::snprintf(buf, sizeof buf, "V3 1 2 3 326.6 0 50\nIM 1 2 3 1e5 400 50 4 %g %g %g %g %g %g %g   # M\n",
                      rs, xls, rr, xlr, xm, H, Tl);
        const auto net = runText(buf, 30000);
        const double sp = currentTrace(*net, "M.speed[pu]").back();
        // equivalent circuit (per unit): Te(s) = |Vth|^2 (rr/s) / ((Rth + rr/s)^2 + (Xth + xlr)^2), load Tl (1-s)^2
        const cd Zm(0, xm), Zs(rs, xls), Vth = Zm / (Zs + Zm);
        const cd Zth = Zs * Zm / (Zs + Zm);
        auto Te = [&](double s) { return std::norm(Vth) * (rr / s) / (std::pow(Zth.real() + rr / s, 2) + std::pow(Zth.imag() + xlr, 2)); };
        double lo = 1e-5, hi = 0.2;   // Te - Tl (1-s)^2 changes sign between
        for (int k = 0; k < 100; ++k) {
            const double m = 0.5 * (lo + hi);
            (Te(m) - Tl * (1 - m) * (1 - m) > 0 ? hi : lo) = m;
        }
        const double s_ec = 0.5 * (lo + hi);
        check(std::fabs((1 - sp) - s_ec) < 0.02 * s_ec + 2e-4, "induction motor: steady-state slip matches the equivalent circuit (slip sim / ec)", (1 - sp) / s_ec);
        const double Ilr = std::abs(1.0 / (Zs + Zm * cd(rr, xlr) / (Zm + cd(rr, xlr))));   // locked-rotor current, pu
        const double Ib = std::sqrt(2.0) * S / (std::sqrt(3.0) * V);
        const auto ia = currentTrace(*net, "M.a");
        const std::vector<double> t6(net->time().begin(), net->time().begin() + 600), i6(ia.begin(), ia.begin() + 600);
        const double i0 = std::abs(phasorOf(t6, i6, fr)) / Ib;   // third 50 Hz cycle after switching on
        check(relErr(i0, Ilr) < 0.15, "induction motor: starting current close to the locked-rotor current (pu, sim)", i0);
    }
}

// ---------------------------------------------------------------- chatter removal and switching interpolation

// Run a netlist with the main program's time loop (Network::runWithEvents when chatter removal
// or interpolation is on, the plain loop otherwise)
std::unique_ptr<Network> runMain(const std::string& text) {
    Circuit ckt = parseNetlist(text, "test");   // sets config from .tran / .options
    const int steps = static_cast<int>(std::lround((config::stpTime - config::srtTime) / config::Dt));
    if (!(config::chatter || config::interp)) return runNetwork(std::move(ckt), steps);
    std::ostringstream sink;
    auto* old = std::cout.rdbuf(sink.rdbuf());
    auto net = std::make_unique<Network>(std::move(ckt.elements), ckt.userNodes);
    net->updateSwitches(config::srtTime);
    net->refactor();
    net->runWithEvents(50);
    std::cout.rdbuf(old);
    return net;
}

void chatterTests() {
    const std::string plain = ".options chatter=off interp=off\n";

    // Watson & Arrillaga, Fig 9.3: opening a switch on an inductor. The trapezoidal rule leaves a
    // sustained step-by-step oscillation of the inductor voltage; chatter removal takes it out.
    const std::string rlOpen = ".tran 0.0001 0.005\nV 0 1 10 0 0\nSW 1 2 closed 0.003\nR 2 3 1\nL 3 0 0.001\n";
    {
        const auto p = runMain(plain + rlOpen), c = runMain(rlOpen);
        const auto vp = nodeTrace(*p, 3), vc = nodeTrace(*c, 3);
        check(p->time().size() == 50 && c->time().size() == 50, "chatter: same output time grid with and without", 0.0);
        check(peakAbs(vp, 35) > 150.0, "chatter: plain trapezoidal rule chatters after opening (V)", peakAbs(vp, 35));
        check(peakAbs(vc, 31) < 1e-3 * peakAbs(vp, 35), "chatter: removed after opening the switch (V)", peakAbs(vc, 31));
    }
    // Fig 9.4: diode-fed RL load, chatter when the diode turns off
    const std::string dRL = ".tran 0.00005 0.03\nV 0 1 100 0 50\nD 1 2 # D\nR 2 3 1\nL 3 0 0.1\n";
    {
        const auto p = runMain(plain + dRL), c = runMain(dRL);
        const auto ip = currentTrace(*p, "D"), vp = nodeTrace(*p, 3), vc = nodeTrace(*c, 3);
        std::size_t k = 1;
        while (k < ip.size() && !(ip[k - 1] > 1e-6 && ip[k] <= 1e-6)) ++k;   // first turn-off
        check(k < ip.size() && peakAbs(vp, k + 2, k + 40) > 10.0, "chatter: diode-fed RL chatters at turn-off (plain, V)",
              peakAbs(vp, k + 2, k + 40));
        check(peakAbs(vc, k + 3, k + 40) < 1e-3 * peakAbs(vp, k + 2, k + 40), "chatter: none after the diode turns off (V)",
              peakAbs(vc, k + 3, k + 40));
    }
    // Fig 9.5: forced commutation. The switch opens between time points; the free-wheeling diode
    // takes over the inductor current at the same instant, so there is no voltage spike.
    {
        const auto c = runMain(".tran 0.00005 0.03\nV 0 1 10 0 0\nSW 1 2 closed 0.020025\nD 0 2\nR 2 3 1\nL 3 0 0.1 # L\n");
        const auto v2 = nodeTrace(*c, 2), iL = currentTrace(*c, "L");
        check(peakAbs(v2, 401) < 0.05, "interpolation: no spike when the free-wheeling diode takes over (V)", peakAbs(v2, 401));
        const double i0 = iL[399], i1 = iL[401];
        check(std::fabs(i1 - i0) < 0.01 * i0, "interpolation: inductor current continues through the diode", std::fabs(i1 - i0));
    }
    // A switch closing between time points (t = 1.02 ms, dt = 0.1 ms) on an RL circuit: with
    // interpolation the current follows i = 10 (1 - exp(-(t - ts)/tau)) closely; without, the switch
    // acts at the next time point (0.08 ms late).
    {
        const std::string rl = ".tran 0.0001 0.005\nV 0 1 10 0 0\nSW 1 2 open 0.00102\nR 2 3 1\nL 3 0 0.001 # L\n";
        auto err = [&](const Network& n) {
            const auto i = currentTrace(n, "L");
            double e = 0;
            for (std::size_t k = 0; k < i.size(); ++k) {
                const double t = n.time()[k], ex = t > 0.00102 ? 10.0 * (1.0 - std::exp(-(t - 0.00102) / 1e-3)) : 0.0;
                e = std::max(e, std::fabs(i[k] - ex));
            }
            return e;
        };
        const auto p = runMain(plain + rl), c = runMain(rl);
        check(err(*c) < 0.4 * err(*p), "interpolation: switching between time points (error A, interp vs plain)", err(*c));
        check(err(*c) < 0.1, "interpolation: error as small as with the switch on a time point (A)", err(*c));
    }
    // A breaker opening at current zero: the current is interrupted at its zero crossing
    {
        const auto c = runMain(".tran 0.00005 0.06\nV 0 1 100 0 50\nSW 1 2 closed 0.0201 czero # BRK\nR 2 3 1\nL 3 0 0.01\n");
        const auto i = currentTrace(*c, "BRK");
        const std::size_t k = static_cast<std::size_t>(std::lround(0.0201 / 0.00005));
        double before = 0.0;
        for (std::size_t q = k - 40; q < k; ++q) before = std::max(before, std::fabs(i[q]));
        check(before > 5.0 && peakAbs(i, k + 300) < 1e-3, "interpolation: breaker opens at a current zero (A after)",
              peakAbs(i, k + 300));
        check(peakAbs(nodeTrace(*c, 3), k + 300) < 0.01, "chatter: none after the breaker opens (V)",
              peakAbs(nodeTrace(*c, 3), k + 300));
    }
    config::chatter = config::defaultChatter;
    config::interp = config::defaultInterp;
}

}  // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);   // unbuffered, so output survives a crash
    const double R = 2.0, L = 10e-3, C = 2e-3, tol = 1e-3;

    {
        auto mk = [&] { Circuit k; k.addSource(0, 1, Vm, a_deg, f); k.addR(1, 2, R); k.addL(0, 2, L); return k; };
        const double e = maxErr(simulate(mk(), 2), [&](double t) { return rl_vL(t, R, L); });
        check(e < tol, "RL, ideal source, L ground->node (max err pu)", e);
    }
    {
        auto mk = [&] { Circuit k; k.addSource(0, 1, Vm, a_deg, f, R); k.addL(1, 0, L); return k; };
        const double e = maxErr(simulate(mk(), 1), [&](double t) { return rl_vL(t, R, L); });
        check(e < tol, "RL, Norton source, L node->ground (max err pu)", e);
    }
    {
        auto mk = [&] { Circuit k; k.addSource(0, 2, Vm, a_deg, f); k.addR(2, 1, R); k.addC(0, 1, C); return k; };
        const double e = maxErr(simulate(mk(), 1), [&](double t) { return rc_vC(t, R, C); });
        check(e < tol, "RC, ideal source on node 2, C ground->node (max err pu)", e);
    }
    {
        auto mk = [&] { Circuit k; k.addSource(1, 0, Vm, a_deg, f); k.addR(1, 2, R); k.addC(2, 0, C); return k; };
        const double e = maxErr(simulate(mk(), 2), [&](double t) { return -rc_vC(t, R, C); });
        check(e < tol, "RC, ideal source reversed (1->0) (max err pu)", e);
    }
    {
        auto mk = [&] { Circuit k; k.addSource(0, 1, Vm, a_deg, f, R); k.addC(0, 1, C); return k; };
        const double e = maxErr(simulate(mk(), 1), [&](double t) { return rc_vC(t, R, C); });
        check(e < tol, "RC, Norton source (max err pu)", e);
    }
    {
        // Switch closing at 20 ms into a series RL circuit
        const double Ron = 0.01, tsw = 0.02, Rt = R + Ron;
        auto mk = [&] {
            Circuit k;
            k.addSource(0, 1, Vm, a_deg, f); k.addSwitch(1, 2, false, {tsw}, Ron);
            k.addR(2, 3, R); k.addL(3, 0, L);
            return k;
        };
        const Run r = simulate(mk(), 3);
        auto exact = [&](double t) {
            const double tp = t - tsw;
            if (tp < 0) return 0.0;
            const cd I = Vm * std::polar(1.0, w * tsw + a_deg * PI / 180.0) / (Rt + j * w * L);
            const double i = std::real(I * std::exp(j * w * tp)) - std::real(I) * std::exp(-Rt * tp / L);
            return E(t) - Rt * i;
        };
        const double e = maxErr(r, exact, tsw + 0.5 * config::Dt);
        check(e < tol, "switch closing into RL at 20 ms (max err pu)", e);
        double pre = 0.0;
        for (std::size_t i = 0; i < r.t.size(); ++i)
            if (r.t[i] < tsw - 0.5 * config::Dt) pre = std::max(pre, std::fabs(r.v[i]) / Vm);
        check(pre < 1e-3, "node isolated while switch open (max |V| pu)", pre);
    }
    {
        Switch sw(1, 2, true, {0.08, 0.03});
        const bool ok = sw.stateAt(0.0) && sw.stateAt(0.0299) && !sw.stateAt(0.03) && !sw.stateAt(0.05) &&
                        sw.stateAt(0.08) && sw.stateAt(0.09);
        check(ok, "switch toggle sequence: open @30 ms, reclose @80 ms", 0.0);
        Fault fl(1, 0, 0.02, 0.04);
        const bool okf = !fl.stateAt(0.0199) && fl.stateAt(0.02) && fl.stateAt(0.0399) && !fl.stateAt(0.04);
        check(okf, "fault applied @20 ms, cleared @40 ms", 0.0);
    }
    {
        // Norton source with internal R == ideal source + separate series R
        auto net = [](Circuit k) {
            k.addL(3, 2, 1e-3); k.addC(2, 1, 1e-3); k.addR(0, 2, 1e-3); k.addL(0, 1, 1e-3);
            return k;
        };
        Circuit a; a.addSource(0, 3, 100, 45, 100, 1.0);
        Circuit b; b.addSource(0, 4, 100, 45, 100); b.addR(4, 3, 1.0);
        const Run ra = simulate(net(std::move(a)), 2), rb = simulate(net(std::move(b)), 2);
        double d = 0.0;
        for (std::size_t i = 0; i < ra.v.size(); ++i) d = std::max(d, std::fabs(ra.v[i] - rb.v[i]));
        check(d < 1e-9, "Norton source == ideal source + series R (max diff V)", d);
    }

    {
        // Sparse LDL^T against dense LU on conductance-like matrices: a meshed grid of
        // branches with a range of conductances (1e-6 .. 1e3 S) and some branches to ground
        auto build = [](std::size_t n, unsigned seed) {
            SparseMatrix A(n, n);
            unsigned s = seed;
            auto rnd = [&s] { s = s * 1103515245u + 12345u; return (s >> 8) / double(1u << 24); };
            auto branch = [&](std::size_t p, std::size_t q, double g) {
                A.add(p, p, g); A.add(q, q, g); A.add(p, q, -g); A.add(q, p, -g);
            };
            const std::size_t width = 10;
            for (std::size_t i = 0; i < n; ++i) {
                const double g = std::pow(10.0, -6.0 + 9.0 * rnd());
                if (i + 1 < n && (i + 1) % width) branch(i, i + 1, g);
                if (i + width < n) branch(i, i + width, std::pow(10.0, -2.0 + 4.0 * rnd()));
                if (i % 7 == 0) A.add(i, i, 0.5 + rnd());   // to ground
            }
            A.finalize();
            return A;
        };
        double worst = 0.0;
        for (std::size_t n : {1u, 2u, 17u, 120u, 400u}) {
            const SparseMatrix A = build(n, 7u + unsigned(n));
            Vector b(n);
            for (std::size_t i = 0; i < n; ++i) b[i] = std::sin(1.0 + double(i));
            SparseLDL ldl; ldl.factor(A);
            LUFactor lu; lu.factor(A.toDense());
            const Vector xs = ldl.solve(b), xd = lu.solve(b);
            double num = 0.0, den = 0.0;
            for (std::size_t i = 0; i < n; ++i) { num = std::max(num, std::fabs(xs[i] - xd[i])); den = std::max(den, std::fabs(xd[i])); }
            worst = std::max(worst, num / den);
        }
        check(worst < 1e-9, "sparse LDL^T matches dense LU (max relative diff)", worst);

        // The 1013-node network's matrix: sparse vs dense, and how little fill-in the ordering gives
        const char* paths[] = {"circuits/grid1000.net", "../circuits/grid1000.net"};
        for (const char* p : paths) {
            std::ifstream probe(p);
            if (!probe) continue;
            std::ostringstream quiet;
            auto* old = std::cout.rdbuf(quiet.rdbuf());
            Circuit gc = loadNetlist(p);
            Network net(std::move(gc.elements), gc.userNodes);
            net.updateSwitches(0.0); net.formG(); net.Gdivide(); net.LU();
            std::cout.rdbuf(old);
            check(net.nnzFactor() < 4 * net.nnzG_UU(), "grid1000: fill-in of L stays small (non-zeros in L)", double(net.nnzFactor()));
            break;
        }

        // Singular matrix (a floating pair of nodes) is reported
        SparseMatrix S(2, 2);
        S.add(0, 0, 1.0); S.add(1, 1, 1.0); S.add(0, 1, -1.0); S.add(1, 0, -1.0); S.finalize();
        bool threw = false;
        try { SparseLDL fac; fac.factor(S); } catch (const std::runtime_error&) { threw = true; }
        check(threw, "sparse LDL^T reports a singular matrix", 0.0);
    }

    try {
        newElementTests();
    } catch (const std::exception& e) {
        std::printf("FAIL  element tests stopped with an error: %s\n", e.what());
        ++failures;
    }
    try {
        chatterTests();
    } catch (const std::exception& e) {
        std::printf("FAIL  chatter tests stopped with an error: %s\n", e.what());
        ++failures;
    }

    std::printf("\n%s\n", failures ? "SOME TESTS FAILED" : "All tests passed");
    return failures ? 1 : 0;
}
