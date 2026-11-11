#include "elements.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "config.h"

namespace {
constexpr double PI = 3.14159265358979323846;
const char* PHASE[] = {"a", "b", "c"};

std::vector<int> concat(std::vector<int> a, const std::vector<int>& b) {
    a.insert(a.end(), b.begin(), b.end());
    return a;
}
std::vector<int> range(int from, int count) {
    std::vector<int> r(count);
    for (int k = 0; k < count; ++k) r[k] = from + k;
    return r;
}
std::string phaseName(const std::string& name, std::size_t k, std::size_t n) {
    return n == 3 ? name + "." + PHASE[k] : n == 1 ? name : name + "." + std::to_string(k + 1);
}
Matrix diag(std::size_t n, double x) {
    Matrix D(n, n);
    for (std::size_t i = 0; i < n; ++i) D(i, i) = x;
    return D;
}
void stamp2(Matrix& g, double G) {
    g(0, 0) += G; g(1, 1) += G; g(0, 1) -= G; g(1, 0) -= G;
}
}  // namespace

// ---------------------------------------------------------------- TwoTerminal

void TwoTerminal::stamp(Matrix& g) const { stamp2(g, 1.0 / Reff); }

void TwoTerminal::history(double t, Vector& h) {
    calcHistory(t);
    h[0] = ihistory;
    h[1] = -ihistory;
}

void TwoTerminal::update(double t, const Vector& v) {
    (void)t;
    Vlast = v[0] - v[1];
    Ilast = Vlast / Reff + ihistory;
}

// ---------------------------------------------------------------- RLC

RLC::RLC(char type, int from, int to, double value) : TwoTerminal(std::string(1, type), from, to), kind_(type) {
    if (!(value > 0.0)) throw std::invalid_argument(std::string(1, type) + ": value must be positive");
    switch (type) {
        case 'R': Reff = value; break;
        case 'L': Reff = 2.0 * value / config::Dt; break;   // 2L/Dt
        case 'C': Reff = config::Dt / (2.0 * value); break; // Dt/2C
        default: throw std::invalid_argument("RLC: type must be R, L or C");
    }
}

void RLC::calcHistory(double t) {
    (void)t;
    switch (kind_) {   // backward-Euler half step: no voltage term
        case 'L': ihistory = config::euler ? Ilast : Ilast + Vlast / Reff; break;
        case 'C': ihistory = config::euler ? -Vlast / Reff : -Ilast - Vlast / Reff; break;
        default: ihistory = 0.0; break;
    }
}

// ---------------------------------------------------------------- Source

Source::Source(int from, int to, double magnitude, double angle_deg, double frequency, double resistance)
    : TwoTerminal("V", from, to), magnitude(magnitude), angle(angle_deg), frequency(frequency), ideal(resistance <= 0.0) {
    if (ideal && from != 0 && to != 0)
        throw std::invalid_argument("An ideal voltage source must have one node connected to ground (node 0)");
    Reff = ideal ? 0.0 : resistance;
}

double Source::E(double t) const { return magnitude * std::cos(2.0 * PI * frequency * t + angle * PI / 180.0); }

void Source::stamp(Matrix& g) const {
    if (!ideal) TwoTerminal::stamp(g);
}

void Source::history(double t, Vector& h) {
    if (ideal) return;
    ihistory = E(t) / Reff;   // Norton equivalent: current source E/R in parallel with R
    h[0] = ihistory;
    h[1] = -ihistory;
}

void Source::update(double t, const Vector& v) {
    if (ideal) Vlast = v[0] - v[1];
    else TwoTerminal::update(t, v);
}

std::vector<int> Source::knownTerminals() const {
    if (!ideal) return {};
    return {nodes[0] == 0 ? 1 : 0};
}

void Source::knownVoltages(double t, Vector& e) const { e[0] = nodes[0] == 0 ? E(t) : -E(t); }

void Source::setKnownCurrents(const Vector& i) { Ilast = nodes[0] == 0 ? i[0] : -i[0]; }

// ---------------------------------------------------------------- Switch

Switch::Switch(int from, int to, bool closed, std::vector<double> t_switch, double Ron, double Roff)
    : Switch("SW", from, to, closed, std::move(t_switch), Ron, Roff) {}

Switch::Switch(std::string type, int from, int to, bool closed, std::vector<double> times, double Ron, double Roff)
    : TwoTerminal(std::move(type), from, to), t_switch(std::move(times)), Ron(Ron), Roff(Roff),
      closed0(closed), closed(closed) {
    std::sort(t_switch.begin(), t_switch.end());
    Reff = closed ? Ron : Roff;
}

// Time accumulates round-off, so a switch time is matched with a small tolerance
double Switch::tol() const { return 1e-6 * config::Dt; }

bool Switch::stateAt(double t) const {
    const auto n = std::count_if(t_switch.begin(), t_switch.end(), [&](double ts) { return t >= ts - tol(); });
    return closed0 != (n % 2 == 1);
}

void Switch::setState(bool c, double t) {
    closed = c;
    Reff = closed ? Ron : Roff;
    pending_open_ = zero_seen_ = false;
    std::printf("t = %.5f s: %s %s (nodes %d-%d)\n", t, label(), describe(closed), nodes[0], nodes[1]);
}

void Switch::commandOpen(double t) {
    pending_open_ = true;
    zero_seen_ = false;
    t_cmd_ = t;
    std::printf("t = %.5f s: %s open command (nodes %d-%d), waiting for current zero\n", t, label(), nodes[0], nodes[1]);
}

bool Switch::timeEvent(double t) {
    const bool want = stateAt(t);
    if (want == closed) {
        if (want) pending_open_ = false;
        return false;
    }
    if (want || !opens_at_current_zero) {   // closing is immediate; so is opening without current-zero detection
        setState(want, t);
        return true;
    }
    if (!pending_open_) commandOpen(t);
    if (!zero_seen_) return false;
    setState(false, t);
    return true;
}

void Switch::update(double t, const Vector& v) {
    const double iprev = Ilast;
    TwoTerminal::update(t, v);
    if (pending_open_ && closed && (iprev * Ilast <= 0.0 || std::fabs(Ilast) <= chop)) zero_seen_ = true;
}

// ---- events (chatter removal / interpolation)
void Switch::beginStep(double a, double b) {
    ia_ = Ilast;
    if (!opens_at_current_zero) return;
    const bool want = stateAt(b);
    if (closed && !want && !pending_open_) {
        double tc = a;   // the command time: the latest switch time up to b
        for (double ts : t_switch)
            if (ts <= b + tol()) tc = std::max(tc, ts);
        commandOpen(tc);
    } else if (want && pending_open_) {
        pending_open_ = false;
    }
}

double Switch::eventAlpha(double a, double dt) const {
    double best = -1.0;
    for (double ts : t_switch)                     // scheduled closing / opening
        if (a + tol() < ts && ts <= a + dt + tol()) {
            const bool want = stateAt(ts);
            if (want != closed && (want || !opens_at_current_zero)) {
                best = std::clamp((ts - a) / dt, 0.0, 1.0);
                break;
            }
        }
    if (pending_open_ && closed) {                 // breaker: current zero (or chopping level)
        const double ia = ia_, ib = Ilast;
        double al = -1.0;
        if (ia == 0.0) al = 0.0;
        else if (ia * ib <= 0.0) al = ia / (ia - ib);
        if (chop > 0.0 && std::fabs(ib) <= chop) {
            const double ac = std::fabs(ia) <= chop ? 0.0 : (std::fabs(ia) - chop) / (std::fabs(ia) - std::fabs(ib));
            al = al < 0.0 ? ac : std::min(al, ac);
        }
        if (al >= 0.0) {
            al = std::clamp(std::max(al, (t_cmd_ - a) / dt), 0.0, 1.0);
            best = best < 0.0 ? al : std::min(best, al);
        }
    }
    return best;
}

bool Switch::applyEvent(double t) {
    const bool want = stateAt(t);
    if (want == closed || (!want && opens_at_current_zero && !pending_open_)) return false;
    setState(want, t);
    return true;
}

// ---------------------------------------------------------------- Fault

namespace {
std::vector<double> faultTimes(double t_apply, std::optional<double> t_clear) {
    if (t_clear) return {t_apply, *t_clear};
    return {t_apply};
}
}  // namespace

Fault::Fault(int from, int to, double t_apply, std::optional<double> t_clear, double Ron, double Roff)
    : Switch("F", from, to, false, faultTimes(t_apply, t_clear), Ron, Roff), t_apply(t_apply), t_clear(t_clear) {}

// ---------------------------------------------------------------- coupled R-L

Matrix sequenceMatrix(double X1, double X0) {
    Matrix M(3, 3);
    const double s = (X0 + 2.0 * X1) / 3.0, m = (X0 - X1) / 3.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) M(i, j) = i == j ? s : m;
    return M;
}

RLBlock::RLBlock(const Matrix& R, const Matrix& L)
    : n(R.rows()), A(inverse(R + (2.0 / config::Dt) * L)), B((2.0 / config::Dt) * L + (-1.0) * R),
      L2((2.0 / config::Dt) * L), vprev(R.rows(), 0.0), iprev(R.rows(), 0.0), ih(R.rows(), 0.0) {}

void RLBlock::stamp(Matrix& g, const std::vector<int>& from, const std::vector<int>& to) const {
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            g(from[i], from[j]) += A(i, j);
            g(to[i], to[j]) += A(i, j);
            g(from[i], to[j]) -= A(i, j);
            g(to[i], from[j]) -= A(i, j);
        }
}

void RLBlock::history(Vector& h, const std::vector<int>& from, const std::vector<int>& to) {
    Vector x(n);
    if (config::euler) {
        x = L2 * iprev;
    } else {
        const Vector Bi = B * iprev;
        for (std::size_t i = 0; i < n; ++i) x[i] = vprev[i] + Bi[i];
    }
    ih = A * x;
    for (std::size_t i = 0; i < n; ++i) {
        h[from[i]] += ih[i];
        h[to[i]] -= ih[i];
    }
}

Vector RLBlock::branchVoltage(const Vector& v, const std::vector<int>& from, const std::vector<int>& to) const {
    Vector vb(n);
    for (std::size_t i = 0; i < n; ++i) vb[i] = v[from[i]] - v[to[i]];
    return vb;
}

void RLBlock::update(const Vector& vb) {
    const Vector Av = A * vb;
    for (std::size_t i = 0; i < n; ++i) iprev[i] = Av[i] + ih[i];
    vprev = vb;
}

CoupledRL::CoupledRL(std::vector<int> from, std::vector<int> to, const Matrix& R, const Matrix& L)
    : Element("RL", concat(from, to)), from_(range(0, int(from.size()))), to_(range(int(from.size()), int(to.size()))),
      rl_(R, L) {
    if (from.size() != to.size() || R.rows() != from.size()) throw std::invalid_argument("CoupledRL: size mismatch");
}

CoupledRL CoupledRL::fromSequence(std::vector<int> from, std::vector<int> to, double R1, double L1, double R0, double L0) {
    return CoupledRL(std::move(from), std::move(to), sequenceMatrix(R1, R0), sequenceMatrix(L1, L0));
}

namespace {
void pushAll(Vector& out, const Vector& x) { out.insert(out.end(), x.begin(), x.end()); }
std::size_t readAll(Vector& x, const Vector& s, std::size_t pos) {
    for (double& v : x) v = s[pos++];
    return pos;
}
}  // namespace

void CoupledRL::saveState(Vector& out) const { pushAll(out, rl_.vprev); pushAll(out, rl_.iprev); }
std::size_t CoupledRL::loadState(const Vector& s, std::size_t pos) { return readAll(rl_.iprev, s, readAll(rl_.vprev, s, pos)); }
void CoupledC::saveState(Vector& out) const { pushAll(out, vprev_); pushAll(out, iprev_); }
std::size_t CoupledC::loadState(const Vector& s, std::size_t pos) { return readAll(iprev_, s, readAll(vprev_, s, pos)); }

std::vector<std::string> CoupledRL::currentNames() const {
    std::vector<std::string> names;
    for (std::size_t k = 0; k < from_.size(); ++k) names.push_back(phaseName(name, k, from_.size()));
    return names;
}

// ---------------------------------------------------------------- coupled shunt C

CoupledC::CoupledC(std::vector<int> nodes, const Matrix& C)
    : Element("C", std::move(nodes)), G_((2.0 / config::Dt) * C), vprev_(C.rows(), 0.0), iprev_(C.rows(), 0.0),
      ih_(C.rows(), 0.0) {}

CoupledC CoupledC::fromSequence(std::vector<int> nodes, double C1, double C0) {
    return CoupledC(std::move(nodes), sequenceMatrix(C1, C0));
}

void CoupledC::stamp(Matrix& g) const {
    for (std::size_t i = 0; i < G_.rows(); ++i)
        for (std::size_t j = 0; j < G_.cols(); ++j) g(i, j) += G_(i, j);
}

void CoupledC::history(double t, Vector& h) {
    (void)t;
    const Vector Gv = G_ * vprev_;
    for (std::size_t i = 0; i < Gv.size(); ++i) {
        ih_[i] = config::euler ? -Gv[i] : -(Gv[i] + iprev_[i]);
        h[i] = ih_[i];
    }
}

void CoupledC::update(double t, const Vector& v) {
    (void)t;
    const Vector Gv = G_ * v;
    for (std::size_t i = 0; i < Gv.size(); ++i) iprev_[i] = Gv[i] + ih_[i];
    vprev_ = v;
}

std::vector<std::string> CoupledC::currentNames() const {
    std::vector<std::string> names;
    for (std::size_t k = 0; k < size(); ++k) names.push_back(phaseName(name, k, size()));
    return names;
}

// ---------------------------------------------------------------- Bergeron line

BergeronLine::BergeronLine(std::vector<int> kEnd, std::vector<int> mEnd, std::vector<double> Zc, std::vector<double> tau,
                           std::vector<double> R)
    : Element("TL", concat(kEnd, mEnd)), M_(kEnd.size()), T_(kEnd.size(), kEnd.size()), tau_(std::move(tau)) {
    if (mEnd.size() != M_ || Zc.size() != M_ || tau_.size() != M_ || R.size() != M_ || (M_ != 1 && M_ != 3))
        throw std::invalid_argument("Bergeron line: 1 or 3 phases, with Zc, tau and R per mode");
    if (M_ == 1) {
        T_(0, 0) = 1.0;
    } else {   // orthonormal Clarke transformation: mode 0 (zero sequence), modes 1 and 2 (aerial)
        const double a = 1.0 / std::sqrt(3.0), b = std::sqrt(2.0 / 3.0), c = 1.0 / std::sqrt(6.0), d = 1.0 / std::sqrt(2.0);
        const double t[3][3] = {{a, a, a}, {b, -c, -c}, {0.0, d, -d}};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) T_(i, j) = t[i][j];
    }
    Matrix Ym(M_, M_);
    for (std::size_t j = 0; j < M_; ++j) {
        if (!(Zc[j] > 0.0)) throw std::invalid_argument("Bergeron line: surge impedance must be positive");
        if (tau_[j] < config::Dt * (1.0 - 1e-9))
            throw std::invalid_argument("Bergeron line: travel time must be at least one time step");
        Z_.push_back(Zc[j] + R[j] / 4.0);
        hk_.push_back((Zc[j] - R[j] / 4.0) / (Zc[j] + R[j] / 4.0));
        Ym(j, j) = 1.0 / Z_[j];
    }
    Yp_ = transpose(T_) * Ym * T_;
    ts_.assign(1, config::srtTime);
    vk_.assign(M_, {0.0}); ik_.assign(M_, {0.0}); vm_.assign(M_, {0.0}); im_.assign(M_, {0.0});
    Ihk_.assign(M_, 0.0); Ihm_.assign(M_, 0.0);
    ikPhase_.assign(M_, 0.0); imPhase_.assign(M_, 0.0);
}

void BergeronLine::stamp(Matrix& g) const {
    for (std::size_t i = 0; i < M_; ++i)
        for (std::size_t j = 0; j < M_; ++j) {
            g(i, j) += Yp_(i, j);
            g(M_ + i, M_ + j) += Yp_(i, j);
        }
}

double BergeronLine::delayed(const std::vector<double>& s, double t, double tau) const {
    const double x = t - tau;
    if (x <= ts_[0]) return x < ts_[0] ? 0.0 : s[0];
    const std::size_t i = static_cast<std::size_t>(std::upper_bound(ts_.begin(), ts_.end(), x) - ts_.begin()) - 1;
    if (i + 1 >= ts_.size()) return s.back();
    const double f = (x - ts_[i]) / (ts_[i + 1] - ts_[i]);
    return s[i] * (1.0 - f) + s[i + 1] * f;
}

void BergeronLine::history(double t, Vector& h) {
    for (std::size_t j = 0; j < M_; ++j) {
        const double Z = Z_[j], hh = hk_[j], tau = tau_[j];
        const double vkd = delayed(vk_[j], t, tau), ikd = delayed(ik_[j], t, tau);
        const double vmd = delayed(vm_[j], t, tau), imd = delayed(im_[j], t, tau);
        Ihk_[j] = 0.5 * (1.0 + hh) * (-vmd / Z - hh * imd) + 0.5 * (1.0 - hh) * (-vkd / Z - hh * ikd);
        Ihm_[j] = 0.5 * (1.0 + hh) * (-vkd / Z - hh * ikd) + 0.5 * (1.0 - hh) * (-vmd / Z - hh * imd);
    }
    const Matrix Tt = transpose(T_);
    const Vector hk = Tt * Ihk_, hm = Tt * Ihm_;
    for (std::size_t p = 0; p < M_; ++p) {
        h[p] += hk[p];
        h[M_ + p] += hm[p];
    }
}

void BergeronLine::update(double t, const Vector& v) {
    Vector vkp(M_), vmp(M_);
    for (std::size_t p = 0; p < M_; ++p) { vkp[p] = v[p]; vmp[p] = v[M_ + p]; }
    const Vector vkm = T_ * vkp, vmm = T_ * vmp;
    Vector ikm(M_), imm(M_);
    for (std::size_t j = 0; j < M_; ++j) {
        ikm[j] = vkm[j] / Z_[j] + Ihk_[j];
        imm[j] = vmm[j] / Z_[j] + Ihm_[j];
        vk_[j].push_back(vkm[j]); ik_[j].push_back(ikm[j]);
        vm_[j].push_back(vmm[j]); im_[j].push_back(imm[j]);
    }
    ts_.push_back(t);
    const Matrix Tt = transpose(T_);
    ikPhase_ = Tt * ikm;
    imPhase_ = Tt * imm;
}

void BergeronLine::saveState(Vector& out) const {
    out.push_back(static_cast<double>(ts_.size()));
    out.push_back(ts_.back());
    for (std::size_t j = 0; j < M_; ++j) {
        out.push_back(vk_[j].back()); out.push_back(ik_[j].back());
        out.push_back(vm_[j].back()); out.push_back(im_[j].back());
    }
    pushAll(out, ikPhase_);
    pushAll(out, imPhase_);
}

std::size_t BergeronLine::loadState(const Vector& s, std::size_t pos) {
    const double c = s[pos];
    const std::size_t k = static_cast<std::size_t>(std::floor(c + 1e-9));
    // keep k samples, then the stored sample becomes the newest: appended if the count is fractional
    // (a state interpolated between two time points), otherwise it replaces sample k-1
    const std::size_t keep = c - static_cast<double>(k) > 1e-9 ? k : k - 1;
    ts_.resize(keep);
    for (std::size_t j = 0; j < M_; ++j)
        for (auto* x : {&vk_[j], &ik_[j], &vm_[j], &im_[j]}) x->resize(keep);
    ts_.push_back(s[pos + 1]);
    for (std::size_t j = 0; j < M_; ++j) {
        const std::size_t q = pos + 2 + 4 * j;
        vk_[j].push_back(s[q]); ik_[j].push_back(s[q + 1]);
        vm_[j].push_back(s[q + 2]); im_[j].push_back(s[q + 3]);
    }
    pos += 2 + 4 * M_;
    pos = readAll(ikPhase_, s, pos);
    return readAll(imPhase_, s, pos);
}

std::vector<std::string> BergeronLine::currentNames() const {
    std::vector<std::string> names;
    for (std::size_t p = 0; p < M_; ++p) names.push_back(phaseName(name + ".k", p, M_));
    for (std::size_t p = 0; p < M_; ++p) names.push_back(phaseName(name + ".m", p, M_));
    return names;
}

void BergeronLine::currents(Vector& out) const {
    for (double i : ikPhase_) out.push_back(i);
    for (double i : imPhase_) out.push_back(i);
}

// ---------------------------------------------------------------- transformer core

TransformerCore::TransformerCore(int p, int q, int s, int t, double n, double R, double L)
    : Element("XF", {p, q, s, t}), a_{1.0, -1.0, -n, n}, n_(n) {
    if (!(n > 0.0)) throw std::invalid_argument("Transformer: ratio must be positive");
    if (!(R + 2.0 * L / config::Dt > 0.0)) throw std::invalid_argument("Transformer: leakage R and L cannot both be zero");
    A_ = 1.0 / (R + 2.0 * L / config::Dt);
    B_ = 2.0 * L / config::Dt - R;
    L2_ = 2.0 * L / config::Dt;
}

void TransformerCore::stamp(Matrix& g) const {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) g(i, j) += A_ * a_[i] * a_[j];
}

void TransformerCore::history(double t, Vector& h) {
    (void)t;
    ih_ = config::euler ? A_ * L2_ * i_ : A_ * (vprev_ + B_ * i_);
    for (int j = 0; j < 4; ++j) h[j] += a_[j] * ih_;
}

void TransformerCore::update(double t, const Vector& v) {
    (void)t;
    double vb = 0.0;
    for (int j = 0; j < 4; ++j) vb += a_[j] * v[j];
    i_ = A_ * vb + ih_;
    vprev_ = vb;
}

// ---------------------------------------------------------------- saturable inductor

SaturableInductor::SaturableInductor(int from, int to, std::vector<double> i, std::vector<double> lambda)
    : Element("LSAT", {from, to}) {
    if (i.empty() || i.size() != lambda.size()) throw std::invalid_argument("LSAT: give pairs of current and flux");
    ip_.push_back(0.0);
    lp_.push_back(0.0);
    for (std::size_t k = 0; k < i.size(); ++k) {
        if (!(i[k] > ip_.back()) || !(lambda[k] > lp_.back()))
            throw std::invalid_argument("LSAT: currents and fluxes must be positive and increasing");
        ip_.push_back(i[k]);
        lp_.push_back(lambda[k]);
    }
}

double SaturableInductor::slopeL(int seg) const {
    const int last = int(ip_.size()) - 1;
    const int k = std::min(seg, last - 1);
    return (lp_[k + 1] - lp_[k]) / (ip_[k + 1] - ip_[k]);
}

void SaturableInductor::segmentOf(double lambda, int& seg, int& sign) const {
    const double a = std::fabs(lambda);
    sign = lambda >= 0.0 ? 1 : -1;
    const int last = int(lp_.size()) - 1;
    seg = last;
    for (int k = 0; k < last; ++k)
        if (a <= lp_[k + 1]) { seg = k; break; }
}

void SaturableInductor::stamp(Matrix& g) const { stamp2(g, config::Dt / (2.0 * slopeL(seg_))); }

// flux before this step's own voltage: lambda_prev + dt/2 v_prev (trapezoidal), lambda_prev (Euler)
double SaturableInductor::flux0() const { return config::euler ? lambda_ : lambda_ + 0.5 * config::Dt * vprev_; }

void SaturableInductor::history(double t, Vector& h) {
    (void)t;
    // segment: i = s i0 + (lambda - s lambda0) / L, with lambda = flux0 + dt/2 v
    const double L = slopeL(seg_), s = sign_;
    ih_ = s * ip_[seg_] + (flux0() - s * lp_[seg_]) / L;
    h[0] = ih_;
    h[1] = -ih_;
}

bool SaturableInductor::check(double t, const Vector& v) {
    (void)t;
    const double lam = flux0() + 0.5 * config::Dt * (v[0] - v[1]);
    int seg = 0, sign = 1;
    segmentOf(lam, seg, sign);
    if (seg != seg_ || (seg != 0 && sign != sign_)) {
        seg_ = seg;
        sign_ = sign;
        return true;
    }
    return false;
}

void SaturableInductor::update(double t, const Vector& v) {
    (void)t;
    const double vb = v[0] - v[1];
    lambda_ = flux0() + 0.5 * config::Dt * vb;
    i_ = config::Dt / (2.0 * slopeL(seg_)) * vb + ih_;
    vprev_ = vb;
}

// ---------------------------------------------------------------- surge arrester

SurgeArrester::SurgeArrester(int from, int to, double Vref, double Iref, double q) : Element("MOV", {from, to}) {
    if (!(Vref > 0.0) || !(Iref > 0.0) || !(q >= 1.0)) throw std::invalid_argument("MOV: need Vref > 0, Iref > 0, q >= 1");
    // Break points from 0.5 Vref to 2.5 Vref, spaced so that the current grows by 30 % per
    // segment (about 1 % linearisation error), plus a first segment from the origin
    vp_.push_back(0.0);
    ip_.push_back(0.0);
    const double ratio = std::pow(1.3, 1.0 / q);
    for (double x = 0.5; x < 2.5 * ratio; x *= ratio) {
        vp_.push_back(x * Vref);
        ip_.push_back(Iref * std::pow(x, q));
    }
}

double SurgeArrester::g(int seg) const {
    const int last = int(vp_.size()) - 1;
    const int k = std::min(seg, last - 1);
    return (ip_[k + 1] - ip_[k]) / (vp_[k + 1] - vp_[k]);
}

double SurgeArrester::i0(int seg, int sign) const { return sign * (ip_[seg] - g(seg) * vp_[seg]); }

void SurgeArrester::segmentOf(double v, int& seg, int& sign) const {
    const double a = std::fabs(v);
    sign = v >= 0.0 ? 1 : -1;
    const int last = int(vp_.size()) - 1;
    seg = last;
    for (int k = 0; k < last; ++k)
        if (a <= vp_[k + 1]) { seg = k; break; }
}

void SurgeArrester::stamp(Matrix& gm) const { stamp2(gm, g(seg_)); }

void SurgeArrester::history(double t, Vector& h) {
    (void)t;
    const double c = i0(seg_, sign_);
    h[0] = c;
    h[1] = -c;
}

bool SurgeArrester::check(double t, const Vector& v) {
    (void)t;
    int seg = 0, sign = 1;
    segmentOf(v[0] - v[1], seg, sign);
    if (seg != seg_ || (seg != 0 && sign != sign_)) {
        seg_ = seg;
        sign_ = sign;
        return true;
    }
    return false;
}

void SurgeArrester::update(double t, const Vector& v) {
    (void)t;
    i_ = g(seg_) * (v[0] - v[1]) + i0(seg_, sign_);
}

// ---------------------------------------------------------------- diode and thyristor

Diode::Diode(int anode, int cathode, double Ron, double Roff, double Vf)
    : Element("D", {anode, cathode}), Ron_(Ron), Roff_(Roff), Vf_(Vf) {
    if (!(Ron > 0.0) || !(Roff > Ron)) throw std::invalid_argument("Diode: need 0 < Ron < Roff");
}

void Diode::stamp(Matrix& g) const { stamp2(g, on_ ? 1.0 / Ron_ : 1.0 / Roff_); }

void Diode::history(double t, Vector& h) {
    (void)t;
    const double c = on_ ? -Vf_ / Ron_ : 0.0;   // on: i = (v - Vf) / Ron
    h[0] = c;
    h[1] = -c;
}

bool Diode::check(double t, const Vector& v) {
    const double vb = v[0] - v[1];
    if (!on_) {
        if (vb > Vf_ && gate(t)) { on_ = true; return true; }
        return false;
    }
    if ((vb - Vf_) / Ron_ < 0.0) { on_ = false; return true; }   // current reversed
    return false;
}

void Diode::update(double t, const Vector& v) {
    (void)t;
    const double vb = v[0] - v[1];
    i_ = on_ ? (vb - Vf_) / Ron_ : vb / Roff_;
    v_ = vb;
}

// turn on where the voltage crosses Vf (with the gate fired), off where the current crosses zero
double Diode::eventAlpha(double a, double dt) const {
    const double b = a + dt;
    if (!on_) {
        const double vb = v_;
        if (vb <= Vf_) return -1.0;
        const double va = va_;
        double al = va >= Vf_ ? 0.0 : (Vf_ - va) / (vb - va);
        const double t = a + al * dt;
        if (!gate(t)) {
            if (!gate(b)) return -1.0;
            double lo = t, hi = b;                 // the gate fires during the step: find when
            for (int k = 0; k < 50; ++k) {
                const double mid = 0.5 * (lo + hi);
                if (gate(mid)) hi = mid; else lo = mid;
            }
            al = (hi - a) / dt;
        }
        return std::clamp(al, 0.0, 1.0);
    }
    const double ib = i_;
    if (ib >= 0.0) return -1.0;
    const double ia = ia_;
    if (ia <= 0.0) return 0.0;
    return std::clamp(ia / (ia - ib), 0.0, 1.0);
}

Thyristor::Thyristor(int anode, int cathode, double alpha_deg, double frequency, double ref_deg, double width_deg,
                     double Ron, double Roff)
    : Diode(anode, cathode, Ron, Roff, 0.0), alpha_(alpha_deg), f_(frequency), ref_(ref_deg), width_(width_deg) {
    type = "THY";
}

bool Thyristor::gate(double t) const {
    const double theta = std::fmod(std::fmod(360.0 * f_ * t + ref_, 360.0) + 360.0, 360.0);
    const double d = std::fmod(theta - alpha_ + 720.0, 360.0);
    return d <= width_;
}

// ---------------------------------------------------------------- synchronous machine

SynchronousMachine::SynchronousMachine(std::vector<int> terminals, std::vector<int> internal, double S, double V, double f,
                                       double H, double xd2, double ra, double E, double delta0_deg, double Pm, double D)
    : Element("SM", concat(terminals, internal)), S_(S), ws_(2.0 * PI * f), Vpk_(V * std::sqrt(2.0 / 3.0)),
      H_(H), E_(E), Pm_(Pm), D_(D), delta_(delta0_deg * PI / 180.0) {
    if (!(S > 0.0) || !(V > 0.0) || !(H > 0.0) || !(xd2 > 0.0)) throw std::invalid_argument("SM: need S, V, H, xd'' > 0");
    const double Zb = V * V / S;
    rl_ = RLBlock(diag(3, ra * Zb), diag(3, xd2 * Zb / ws_));
}

void SynchronousMachine::knownVoltages(double t, Vector& e) const {
    const double th = ws_ * t + delta_;
    for (int k = 0; k < 3; ++k) e[k] = E_ * Vpk_ * std::cos(th - 2.0 * PI * k / 3.0);
}

void SynchronousMachine::update(double t, const Vector& v) {
    rl_.update(rl_.branchVoltage(v, xpos_, tpos_));   // current from the EMF to the terminal (generated)
    Vector e(3);
    knownVoltages(t, e);
    Pe_ = (e[0] * rl_.iprev[0] + e[1] * rl_.iprev[1] + e[2] * rl_.iprev[2]) / S_;
    w_ += config::step * (Pm_ - Pe_ - D_ * (w_ - 1.0)) / (2.0 * H_);   // swing equation (per unit)
    delta_ += config::step * ws_ * (w_ - 1.0);
}

void SynchronousMachine::saveState(Vector& out) const {
    pushAll(out, rl_.vprev); pushAll(out, rl_.iprev);
    out.push_back(w_); out.push_back(delta_); out.push_back(Pe_);
}

std::size_t SynchronousMachine::loadState(const Vector& s, std::size_t pos) {
    pos = readAll(rl_.iprev, s, readAll(rl_.vprev, s, pos));
    w_ = s[pos]; delta_ = s[pos + 1]; Pe_ = s[pos + 2];
    return pos + 3;
}

std::vector<std::string> SynchronousMachine::currentNames() const {
    return {name + ".a", name + ".b", name + ".c", name + ".speed[pu]", name + ".delta[deg]", name + ".Pe[pu]"};
}

void SynchronousMachine::currents(Vector& out) const {
    for (double i : rl_.iprev) out.push_back(i);
    out.push_back(w_);
    out.push_back(delta_ * 180.0 / PI);
    out.push_back(Pe_);
}

// ---------------------------------------------------------------- induction machine

InductionMachine::InductionMachine(std::vector<int> terminals, std::vector<int> internal, double S, double V, double f,
                                   int poles, double rs, double xls, double rr, double xlr, double xm, double H,
                                   double Tload, double speed0)
    : Element("IM", concat(terminals, internal)), Tload_(Tload) {
    if (!(S > 0.0) || !(V > 0.0) || poles < 2 || poles % 2 || !(xm > 0.0) || !(H > 0.0))
        throw std::invalid_argument("IM: need S, V, H, xm > 0 and an even number of poles");
    const double ws = 2.0 * PI * f, Zb = V * V / S;
    P2_ = poles / 2.0;
    wms_ = ws / P2_;
    const double Lls = xls * Zb / ws, Llr = xlr * Zb / ws;
    Lm_ = xm * Zb / ws;
    Lr_ = Llr + Lm_;
    Rr_ = rr * Zb;
    const double Lp = Lls + Lm_ - Lm_ * Lm_ / Lr_;   // transient inductance L'
    rl_ = RLBlock(diag(3, rs * Zb), diag(3, Lp));
    J_ = 2.0 * H * S / (wms_ * wms_);
    Tb_ = S / wms_;
    wm_ = speed0 * wms_;
}

void InductionMachine::knownVoltages(double t, Vector& e) const {
    (void)t;
    const cplx a2 = std::polar(1.0, -2.0 * PI / 3.0), a1 = std::polar(1.0, 2.0 * PI / 3.0);
    e[0] = eprime_.real();
    e[1] = (eprime_ * a2).real();
    e[2] = (eprime_ * a1).real();
}

void InductionMachine::update(double t, const Vector& v) {
    (void)t;
    rl_.update(rl_.branchVoltage(v, tpos_, xpos_));   // current from the terminal into the machine
    const cplx a = std::polar(1.0, 2.0 * PI / 3.0);
    const cplx is = (2.0 / 3.0) * (rl_.iprev[0] + a * rl_.iprev[1] + a * a * rl_.iprev[2]);
    const double h = config::step;
    const cplx A = cplx(-Rr_ / Lr_, P2_ * wm_), B = Rr_ * Lm_ / Lr_;
    const cplx psi = (psir_ * (1.0 + A * h / 2.0) + B * h / 2.0 * (is + isprev_)) / (1.0 - A * h / 2.0);
    Te_ = 1.5 * P2_ * (Lm_ / Lr_) * (psi.real() * is.imag() - psi.imag() * is.real());
    const double Tl = Tload_ * Tb_ * (wm_ / wms_) * (wm_ / wms_);
    wm_ += h * (Te_ - Tl) / J_;
    const cplx A2 = cplx(-Rr_ / Lr_, P2_ * wm_);
    eprime_ = (Lm_ / Lr_) * (A2 * psi + B * is);
    psir_ = psi;
    isprev_ = is;
}

void InductionMachine::saveState(Vector& out) const {
    pushAll(out, rl_.vprev); pushAll(out, rl_.iprev);
    out.push_back(wm_);
    for (const cplx& c : {psir_, isprev_, eprime_}) { out.push_back(c.real()); out.push_back(c.imag()); }
    out.push_back(Te_);
}

std::size_t InductionMachine::loadState(const Vector& s, std::size_t pos) {
    pos = readAll(rl_.iprev, s, readAll(rl_.vprev, s, pos));
    wm_ = s[pos++];
    for (cplx* c : {&psir_, &isprev_, &eprime_}) { *c = cplx(s[pos], s[pos + 1]); pos += 2; }
    Te_ = s[pos];
    return pos + 1;
}

std::vector<std::string> InductionMachine::currentNames() const {
    return {name + ".a", name + ".b", name + ".c", name + ".speed[pu]", name + ".torque[pu]"};
}

void InductionMachine::currents(Vector& out) const {
    for (double i : rl_.iprev) out.push_back(i);
    out.push_back(wm_ / wms_);
    out.push_back(Te_ / Tb_);
}
