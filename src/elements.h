#pragma once

// Circuit elements.
//
// Every element connects k terminals (node numbers, 0 = ground). For each time
// step it behaves, as seen from the network, like
//
//     i = g v + h
//
// where v are the terminal voltages, i the currents flowing from each terminal
// node into the element, g a k x k conductance matrix (stamp) and h a vector of
// history currents (history). This covers single branches (R, L, C), coupled
// branches, lines, transformers and machines alike.
//
// Per time step the network calls, in this order:
//   timeEvent(t)  scheduled state changes (switches)         -> true: g changed
//   history(t, h) history currents for the step
//   ... network solution ...
//   check(t, v)   state changes that depend on the solution   -> true: g changed,
//                 (diodes, arresters, saturation)                  re-solve the step
//   update(t, v)  commit the solution
//
// Ideal voltage sources (and machines' internal EMFs) fix the voltage of some
// terminals: knownTerminals() / knownVoltages(); the network then tells them the
// current they deliver with setKnownCurrents().
//
// Chatter removal and switching interpolation add: a backward-Euler half step when config::euler
// is set (same conductances, no voltage term in the histories; config::step is the step length),
// saveState/loadState for the continuous state, and for switching elements beginStep / eventAlpha
// (where in the step, 0..1, it should have switched; < 0: not) / applyEvent.

#include <complex>
#include <optional>
#include <string>
#include <vector>

#include "linalg.h"

class Element {
public:
    Element(std::string type, std::vector<int> nodes) : type(std::move(type)), nodes(std::move(nodes)) {}
    virtual ~Element() = default;

    std::string type;          // netlist keyword (R, L, SW, TL3, ...)
    std::string name;          // for recorded currents and messages
    std::vector<int> nodes;    // terminal node numbers (0 = ground)
    std::vector<int> idx;      // matrix index of each terminal (-1 = ground), set by the network

    std::size_t size() const { return nodes.size(); }

    virtual void stamp(Matrix& g) const = 0;              // add the k x k conductance (g is zeroed)
    virtual void history(double t, Vector& h) = 0;        // h is zeroed, size k
    virtual bool timeEvent(double t) { (void)t; return false; }
    virtual bool check(double t, const Vector& v) { (void)t; (void)v; return false; }
    virtual bool nonlinear() const { return false; }   // true if check() can change the state
    virtual void update(double t, const Vector& v) = 0;

    virtual std::vector<int> knownTerminals() const { return {}; }      // positions in nodes[]
    virtual void knownVoltages(double t, Vector& e) const { (void)t; (void)e; }
    virtual void setKnownCurrents(const Vector& i) { (void)i; }         // delivered into the network

    virtual std::vector<std::string> currentNames() const { return {name}; }
    virtual void currents(Vector& out) const = 0;                       // append, same order as names

    // chatter removal / switching interpolation
    virtual void saveState(Vector& out) const { (void)out; }            // append the continuous state
    virtual std::size_t loadState(const Vector& s, std::size_t pos) { (void)s; return pos; }
    virtual bool hasEvents() const { return false; }   // switches, faults, diodes, thyristors
    virtual bool valve() const { return false; }       // diodes and thyristors: switched by events, not check()
    virtual bool valveOn() const { return false; }     // a valve that is conducting
    virtual void beginStep(double a, double b) { (void)a; (void)b; }
    virtual double eventAlpha(double a, double dt) const { (void)a; (void)dt; return -1.0; }
    virtual bool applyEvent(double t) { (void)t; return false; }
};

// ---------------------------------------------------------------- two-terminal branches
// i (from nodes[0] to nodes[1]) = V / Reff + ihistory,  V = v[0] - v[1]

class TwoTerminal : public Element {
public:
    TwoTerminal(std::string type, int from, int to) : Element(std::move(type), {from, to}) {}

    void stamp(Matrix& g) const override;
    void history(double t, Vector& h) override;
    void update(double t, const Vector& v) override;
    void currents(Vector& out) const override { out.push_back(Ilast); }
    void saveState(Vector& out) const override { out.push_back(Ilast); out.push_back(Vlast); }
    std::size_t loadState(const Vector& s, std::size_t pos) override { Ilast = s[pos]; Vlast = s[pos + 1]; return pos + 2; }

    int strnode() const { return nodes[0]; }
    int stpnode() const { return nodes[1]; }

    double Reff = 1.0;
    double Ilast = 0.0;     // branch current at the previous time step
    double Vlast = 0.0;     // branch voltage at the previous time step
    double ihistory = 0.0;

protected:
    virtual void calcHistory(double t) { (void)t; ihistory = 0.0; }
};

// Resistor, inductor or capacitor with trapezoidal-rule companion model
class RLC : public TwoTerminal {
public:
    RLC(char type, int from, int to, double value);   // type 'R', 'L' or 'C'
protected:
    void calcHistory(double t) override;
private:
    char kind_;
};

// Sinusoidal voltage source E(t) = magnitude cos(2 pi f t + angle), raising nodes[1] above nodes[0].
//   resistance > 0 : Norton equivalent;  resistance <= 0 : ideal (one terminal on ground)
class Source : public TwoTerminal {
public:
    Source(int from, int to, double magnitude, double angle_deg, double frequency, double resistance = 0.0);

    double E(double t) const;
    void stamp(Matrix& g) const override;
    void history(double t, Vector& h) override;
    void update(double t, const Vector& v) override;
    std::vector<int> knownTerminals() const override;
    void knownVoltages(double t, Vector& e) const override;
    void setKnownCurrents(const Vector& i) override;

    double magnitude, angle, frequency;
    bool ideal;
};

// Time-controlled switch: Ron closed, Roff open; toggles at each time in t_switch.
// With opens_at_current_zero, an opening command takes effect at the next current zero
// (or as soon as |i| <= chop, to model current chopping).
class Switch : public TwoTerminal {
public:
    Switch(int from, int to, bool closed, std::vector<double> t_switch, double Ron = 0.01, double Roff = 1.0e6);

    bool stateAt(double t) const;
    bool timeEvent(double t) override;
    void update(double t, const Vector& v) override;

    bool hasEvents() const override { return true; }
    void beginStep(double a, double b) override;
    double eventAlpha(double a, double dt) const override;
    bool applyEvent(double t) override;

    std::vector<double> t_switch;
    double Ron, Roff;
    bool closed0, closed;
    bool opens_at_current_zero = false;
    double chop = 0.0;

protected:
    Switch(std::string type, int from, int to, bool closed, std::vector<double> t_switch, double Ron, double Roff);
    double tol() const;                 // time tolerance for switch times
    void setState(bool c, double t);    // switch, and forget any pending open command
    void commandOpen(double t);         // wait for the next current zero
    virtual const char* label() const { return "switch"; }
    virtual const char* describe(bool isClosed) const { return isClosed ? "closed" : "opened"; }
    bool pending_open_ = false, zero_seen_ = false;
    double t_cmd_ = 0.0, ia_ = 0.0;
};

// Fault: closed from t_apply until t_clear (no t_clear = permanent)
class Fault : public Switch {
public:
    Fault(int from, int to, double t_apply, std::optional<double> t_clear, double Ron = 0.01, double Roff = 1.0e6);
    double t_apply;
    std::optional<double> t_clear;
protected:
    const char* label() const override { return "fault"; }
    const char* describe(bool isClosed) const override { return isClosed ? "applied" : "cleared"; }
};

// ---------------------------------------------------------------- coupled R-L
// n coupled series branches (from[k] -> to[k]) with resistance and inductance matrices.
// Trapezoidal rule:  i = A v + A [ v_prev + (2L/dt - R) i_prev ],  A = (R + 2L/dt)^-1

struct RLBlock {
    RLBlock() = default;
    RLBlock(const Matrix& R, const Matrix& L);
    std::size_t n = 0;
    Matrix A, B, L2;                  // A = (R + 2L/dt)^-1,  B = 2L/dt - R,  L2 = 2L/dt
    Vector vprev, iprev, ih;
    void stamp(Matrix& g, const std::vector<int>& from, const std::vector<int>& to) const;
    void history(Vector& h, const std::vector<int>& from, const std::vector<int>& to);
    Vector branchVoltage(const Vector& v, const std::vector<int>& from, const std::vector<int>& to) const;
    void update(const Vector& vb);    // commit: i = A vb + ih
};

class CoupledRL : public Element {
public:
    // nodes = from (n) then to (n)
    CoupledRL(std::vector<int> from, std::vector<int> to, const Matrix& R, const Matrix& L);
    // three-phase from sequence data: R1, L1 (positive) and R0, L0 (zero sequence)
    static CoupledRL fromSequence(std::vector<int> from, std::vector<int> to, double R1, double L1, double R0, double L0);

    void stamp(Matrix& g) const override { rl_.stamp(g, from_, to_); }
    void history(double t, Vector& h) override { (void)t; rl_.history(h, from_, to_); }
    void update(double t, const Vector& v) override { (void)t; rl_.update(rl_.branchVoltage(v, from_, to_)); }
    std::vector<std::string> currentNames() const override;
    void currents(Vector& out) const override { for (double i : rl_.iprev) out.push_back(i); }
    void saveState(Vector& out) const override;
    std::size_t loadState(const Vector& s, std::size_t pos) override;

private:
    std::vector<int> from_, to_;      // positions in nodes[]
    RLBlock rl_;
};

// n coupled shunt capacitances to ground (capacitance matrix C)
class CoupledC : public Element {
public:
    CoupledC(std::vector<int> nodes, const Matrix& C);
    static CoupledC fromSequence(std::vector<int> nodes, double C1, double C0);

    void stamp(Matrix& g) const override;
    void history(double t, Vector& h) override;
    void update(double t, const Vector& v) override;
    std::vector<std::string> currentNames() const override;
    void currents(Vector& out) const override { for (double i : iprev_) out.push_back(i); }
    void saveState(Vector& out) const override;
    std::size_t loadState(const Vector& s, std::size_t pos) override;

private:
    Matrix G_;                        // 2C/dt
    Vector vprev_, iprev_, ih_;
};

// Symmetric components -> phase-domain self and mutual values for a 3x3 matrix
Matrix sequenceMatrix(double X1, double X0);

// ---------------------------------------------------------------- Bergeron line
// Travelling-wave line with M modes (M = 1: single-phase; M = 3: transposed three-phase line,
// Clarke modal transformation). Per mode: surge impedance Zc, travel time tau, total series
// resistance R (lumped R/4, R/2, R/4 as in the EMTP Bergeron model). nodes = k end (M) then m end (M).
class BergeronLine : public Element {
public:
    BergeronLine(std::vector<int> kEnd, std::vector<int> mEnd, std::vector<double> Zc, std::vector<double> tau,
                 std::vector<double> R);

    void stamp(Matrix& g) const override;
    void history(double t, Vector& h) override;
    void update(double t, const Vector& v) override;
    std::vector<std::string> currentNames() const override;
    void currents(Vector& out) const override;
    // state: number of samples, time and values of the newest sample, phase currents. Loading a
    // state between two saved ones (a fractional count) replaces the newest sample.
    void saveState(Vector& out) const override;
    std::size_t loadState(const Vector& s, std::size_t pos) override;

private:
    std::size_t M_;
    Matrix T_;                        // phase -> mode (orthonormal)
    std::vector<double> Z_, tau_, hk_;  // Z' = Zc + R/4,  h = (Zc - R/4) / (Zc + R/4)
    Matrix Yp_;                       // phase-domain conductance at each end: T^T diag(1/Z') T
    // modal histories, one sample per solution point at times ts_ (sample 0 = start time)
    std::vector<double> ts_;
    std::vector<std::vector<double>> vk_, ik_, vm_, im_;
    Vector Ihk_, Ihm_;                // modal history currents for this step
    Vector ikPhase_, imPhase_;        // phase currents into the line at each end (last step)
    double delayed(const std::vector<double>& s, double t, double tau) const;
};

// ---------------------------------------------------------------- transformer
// Ideal transformer (ratio n = V1/V2) with series leakage R, L on the primary side.
// nodes = p, q (primary), s, t (secondary). Branch voltage vb = (vp - vq) - n (vs - vt).
// The magnetising branch is a separate element across p-q.
class TransformerCore : public Element {
public:
    TransformerCore(int p, int q, int s, int t, double n, double R, double L);

    void stamp(Matrix& g) const override;
    void history(double t, Vector& h) override;
    void update(double t, const Vector& v) override;
    std::vector<std::string> currentNames() const override { return {name + ".1", name + ".2"}; }
    void currents(Vector& out) const override { out.push_back(i_); out.push_back(n_ * i_); }
    void saveState(Vector& out) const override { out.push_back(vprev_); out.push_back(i_); }
    std::size_t loadState(const Vector& s, std::size_t pos) override { vprev_ = s[pos]; i_ = s[pos + 1]; return pos + 2; }

private:
    double a_[4];
    double n_, A_, B_, L2_;           // A = 1/(R + 2L/dt),  B = 2L/dt - R,  L2 = 2L/dt
    double vprev_ = 0.0, i_ = 0.0, ih_ = 0.0;
};

// ---------------------------------------------------------------- non-linear elements

// Non-linear inductor with a symmetric piecewise-linear flux-current curve through the origin:
// points (i_k, lambda_k), k = 1..n, increasing; beyond the last point the last slope continues.
class SaturableInductor : public Element {
public:
    SaturableInductor(int from, int to, std::vector<double> i, std::vector<double> lambda);

    void stamp(Matrix& g) const override;
    void history(double t, Vector& h) override;
    bool check(double t, const Vector& v) override;
    bool nonlinear() const override { return true; }
    void update(double t, const Vector& v) override;
    void currents(Vector& out) const override { out.push_back(i_); }
    double flux() const { return lambda_; }
    void saveState(Vector& out) const override { out.push_back(lambda_); out.push_back(vprev_); out.push_back(i_); }
    std::size_t loadState(const Vector& s, std::size_t pos) override {
        lambda_ = s[pos]; vprev_ = s[pos + 1]; i_ = s[pos + 2]; return pos + 3;
    }

private:
    std::vector<double> ip_, lp_;     // break points including (0, 0)
    double flux0() const;             // flux before this step's own voltage
    int seg_ = 0, sign_ = 1;
    double lambda_ = 0.0, vprev_ = 0.0, i_ = 0.0, ih_ = 0.0;
    void segmentOf(double lambda, int& seg, int& sign) const;
    double slopeL(int seg) const;
};

// Metal-oxide surge arrester: i = Iref (|v| / Vref)^q sign(v), as a piecewise-linear resistor
class SurgeArrester : public Element {
public:
    SurgeArrester(int from, int to, double Vref, double Iref = 1000.0, double q = 26.0);

    void stamp(Matrix& g) const override;
    void history(double t, Vector& h) override;
    bool check(double t, const Vector& v) override;
    bool nonlinear() const override { return true; }
    void update(double t, const Vector& v) override;
    void currents(Vector& out) const override { out.push_back(i_); }
    void saveState(Vector& out) const override { out.push_back(i_); }
    std::size_t loadState(const Vector& s, std::size_t pos) override { i_ = s[pos]; return pos + 1; }

private:
    std::vector<double> vp_, ip_;     // break points including (0, 0)
    int seg_ = 0, sign_ = 1;
    double i_ = 0.0;
    void segmentOf(double v, int& seg, int& sign) const;
    double g(int seg) const;
    double i0(int seg, int sign) const;   // i = g v + i0
};

// Diode (anode nodes[0], cathode nodes[1]): Ron (in series with Vf) when on, Roff when off.
// Turns on when v > Vf, off when the current reverses.
class Diode : public Element {
public:
    Diode(int anode, int cathode, double Ron = 1e-3, double Roff = 1e6, double Vf = 0.0);

    void stamp(Matrix& g) const override;
    void history(double t, Vector& h) override;
    bool check(double t, const Vector& v) override;
    bool nonlinear() const override { return true; }
    void update(double t, const Vector& v) override;
    void currents(Vector& out) const override { out.push_back(i_); }
    bool on() const { return on_; }
    bool valveOn() const override { return on_; }
    void saveState(Vector& out) const override { out.push_back(i_); out.push_back(v_); }
    std::size_t loadState(const Vector& s, std::size_t pos) override { i_ = s[pos]; v_ = s[pos + 1]; return pos + 2; }
    bool hasEvents() const override { return true; }
    bool valve() const override { return true; }
    void beginStep(double a, double b) override { (void)a; (void)b; ia_ = i_; va_ = v_; }
    double eventAlpha(double a, double dt) const override;
    bool applyEvent(double t) override { (void)t; on_ = !on_; return true; }

protected:
    virtual bool gate(double t) const { (void)t; return true; }
    double Ron_, Roff_, Vf_;
    bool on_ = false;
    double i_ = 0.0, v_ = 0.0, ia_ = 0.0, va_ = 0.0;
};

// Thyristor: a diode that can only turn on while its gate is fired. Simple control:
// the gate is on while the phase (360 f t + ref) mod 360 lies in [alpha, alpha + width] degrees.
class Thyristor : public Diode {
public:
    Thyristor(int anode, int cathode, double alpha_deg, double frequency, double ref_deg = 0.0,
              double width_deg = 120.0, double Ron = 1e-3, double Roff = 1e6);
protected:
    bool gate(double t) const override;
private:
    double alpha_, f_, ref_, width_;
};

// ---------------------------------------------------------------- machines

// Synchronous machine, classical model: EMF E behind subtransient reactance xd'' (and ra),
// rotor angle and speed from the swing equation 2H dw/dt = Pm - Pe - D (w - 1)  (per unit).
// nodes = a, b, c (terminals), xa, xb, xc (internal EMF nodes, created by the netlist reader).
class SynchronousMachine : public Element {
public:
    SynchronousMachine(std::vector<int> terminals, std::vector<int> internal, double S, double V, double f,
                       double H, double xd2, double ra, double E, double delta0_deg, double Pm, double D);

    void stamp(Matrix& g) const override { rl_.stamp(g, xpos_, tpos_); }
    void history(double t, Vector& h) override { (void)t; rl_.history(h, xpos_, tpos_); }
    void update(double t, const Vector& v) override;
    std::vector<int> knownTerminals() const override { return xpos_; }
    void knownVoltages(double t, Vector& e) const override;
    std::vector<std::string> currentNames() const override;
    void currents(Vector& out) const override;
    void saveState(Vector& out) const override;
    std::size_t loadState(const Vector& s, std::size_t pos) override;

    double speed() const { return w_; }
    double delta() const { return delta_; }

private:
    std::vector<int> tpos_{0, 1, 2}, xpos_{3, 4, 5};
    RLBlock rl_;
    double S_, ws_, Vpk_, H_, E_, Pm_, D_;
    double w_ = 1.0, delta_ = 0.0, Pe_ = 0.0;
};

// Induction machine (motor convention), voltage behind transient inductance in the stator frame:
//   v = Rs i + L' di/dt + e',  e' = (Lm/Lr) d(psi_r)/dt,
//   d(psi_r)/dt = -(Rr/Lr)(psi_r - Lm i) + j w_r psi_r,  J dw_m/dt = Te - T_load (T_load ~ w^2).
// Parameters in per unit on the machine rating. nodes = a, b, c, xa, xb, xc (internal).
class InductionMachine : public Element {
public:
    InductionMachine(std::vector<int> terminals, std::vector<int> internal, double S, double V, double f, int poles,
                     double rs, double xls, double rr, double xlr, double xm, double H, double Tload, double speed0);

    void stamp(Matrix& g) const override { rl_.stamp(g, tpos_, xpos_); }
    void history(double t, Vector& h) override { (void)t; rl_.history(h, tpos_, xpos_); }
    void update(double t, const Vector& v) override;
    std::vector<int> knownTerminals() const override { return xpos_; }
    void knownVoltages(double t, Vector& e) const override;
    std::vector<std::string> currentNames() const override;
    void currents(Vector& out) const override;
    void saveState(Vector& out) const override;
    std::size_t loadState(const Vector& s, std::size_t pos) override;

    double speed() const { return wm_ / wms_; }

private:
    using cplx = std::complex<double>;
    std::vector<int> tpos_{0, 1, 2}, xpos_{3, 4, 5};
    RLBlock rl_;
    double Rr_, Lr_, Lm_, P2_, J_, Tb_, Tload_, wms_;
    cplx psir_{0.0, 0.0}, isprev_{0.0, 0.0}, eprime_{0.0, 0.0};
    double wm_ = 0.0, Te_ = 0.0;
};
