#include "network.h"

#include "config.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

Network::Network(std::vector<std::unique_ptr<Element>> elements, int userNodes) : comp_list_(std::move(elements)) {
    for (const auto& e : comp_list_)
        for (int n : e->nodes) {
            if (n < 0) throw std::invalid_argument("Negative node number in " + e->name);
            numNodes_ = std::max(numNodes_, n);
        }
    userNodes_ = userNodes < 0 ? numNodes_ : std::min(userNodes, numNodes_);

    // Known nodes: voltages fixed by ideal sources and machine EMFs
    std::vector<int> known;
    for (const auto& e : comp_list_)
        for (int p : e->knownTerminals()) {
            const int node = e->nodes[p];
            if (node == 0) throw std::invalid_argument(e->name + ": a known (ideal source) voltage cannot be on ground");
            if (std::find(known.begin(), known.end(), node) != known.end())
                throw std::invalid_argument("Two ideal voltage sources on node " + std::to_string(node));
            known.push_back(node);
        }

    // Internal node ordering: unknown nodes first, then known nodes, so that G can
    // be partitioned into G_UU, G_UK, G_KU, G_KK regardless of how the nodes are numbered.
    std::vector<char> isKnown(numNodes_ + 1, 0);
    for (int n : known) isKnown[n] = 1;
    for (int n = 1; n <= numNodes_; ++n)
        if (!isKnown[n]) order_.push_back(n);
    order_.insert(order_.end(), known.begin(), known.end());

    idx_.assign(numNodes_ + 1, -1);
    for (int i = 0; i < numNodes_; ++i) idx_[order_[i]] = i;
    src_count_ = static_cast<int>(known.size());
    unknownNodes_ = numNodes_ - src_count_;

    for (auto& e : comp_list_) {
        e->idx.clear();
        for (int n : e->nodes) e->idx.push_back(idx_[n]);
        std::vector<int> ki;
        for (int p : e->knownTerminals()) ki.push_back(idx_[e->nodes[p]] - unknownNodes_);
        knownIdx_.push_back(ki);
        hist_.emplace_back(e->size(), 0.0);
        if (e->nonlinear()) nonlinear_.push_back(hist_.size() - 1);
        if (e->hasEvents()) events_.push_back(e.get());
        for (const auto& nm : e->currentNames()) currentNames_.push_back(nm);
    }
    Vn_.assign(numNodes_, 0.0);
    I_H_.assign(numNodes_, 0.0);

    std::cout << "No of ideal sources " << src_count_ << "\n";
    std::cout << "No of unknown nodes " << unknownNodes_ << "\n\n";
}

bool Network::updateSwitches(double t) {
    bool changed = false;
    for (auto& e : comp_list_) changed = e->timeEvent(t) || changed;
    return changed;
}

void Network::formG() {
    G_ = SparseMatrix(numNodes_, numNodes_);
    for (const auto& e : comp_list_) {
        const std::size_t k = e->size();
        Matrix g(k, k);
        e->stamp(g);
        for (std::size_t i = 0; i < k; ++i) {
            if (e->idx[i] < 0) continue;
            for (std::size_t j = 0; j < k; ++j)
                if (e->idx[j] >= 0 && g(i, j) != 0.0) G_.add(e->idx[i], e->idx[j], g(i, j));
        }
    }
    G_.finalize();   // entries stamped at the same position are summed
    if (numNodes_ <= 10) G_.toDense().print(std::cout, "G");   // printing a large G would swamp the output
}

void Network::Gdivide() {
    const std::size_t u = unknownNodes_, k = src_count_;
    G_UU_ = G_.block(0, 0, u, u);
    G_UK_ = G_.block(0, u, u, k);
    G_KU_ = G_.block(u, 0, k, u);
    G_KK_ = G_.block(u, u, k, k);
}

void Network::LU() {
    G_UU_LU_.factor(G_UU_);
    ++factorisations_;
}

void Network::Vdivide() {
    V_U_.assign(Vn_.begin(), Vn_.begin() + unknownNodes_);
    V_K_.assign(Vn_.begin() + unknownNodes_, Vn_.end());
}

void Network::Vkupdate(double t) {
    for (std::size_t n = 0; n < comp_list_.size(); ++n) {
        if (knownIdx_[n].empty()) continue;
        Vector e(knownIdx_[n].size(), 0.0);
        comp_list_[n]->knownVoltages(t, e);
        for (std::size_t j = 0; j < e.size(); ++j) V_K_[knownIdx_[n][j]] = e[j];
    }
}

void Network::calcBrnHistory(double t) {
    for (std::size_t n = 0; n < comp_list_.size(); ++n) {
        std::fill(hist_[n].begin(), hist_[n].end(), 0.0);
        comp_list_[n]->history(t, hist_[n]);
    }
}

void Network::calcinjection() {
    I_H_.assign(numNodes_, 0.0);
    // Each element draws i = g v + h from its terminals, so h is subtracted at each terminal node
    for (std::size_t n = 0; n < comp_list_.size(); ++n) {
        const auto& e = comp_list_[n];
        for (std::size_t j = 0; j < e->size(); ++j)
            if (e->idx[j] >= 0) I_H_[e->idx[j]] -= hist_[n][j];
    }
}

void Network::Idivide() {
    I_U_.assign(I_H_.begin(), I_H_.begin() + unknownNodes_);
    I_K_.assign(I_H_.begin() + unknownNodes_, I_H_.end());
}

void Network::calIUnew() {
    I_U_new_ = I_U_;
    if (src_count_ > 0) {
        const Vector gv = G_UK_ * V_K_;
        for (int i = 0; i < unknownNodes_; ++i) I_U_new_[i] -= gv[i];
    }
}

void Network::calcnewV() { V_U_ = G_UU_LU_.solve(I_U_new_); }

void Network::reconstruct_V() {
    std::copy(V_U_.begin(), V_U_.end(), Vn_.begin());
    std::copy(V_K_.begin(), V_K_.end(), Vn_.begin() + unknownNodes_);
}

const Vector& Network::terminalVoltages(const Element& e) {
    vbuf_.resize(e.size());
    for (std::size_t j = 0; j < e.size(); ++j) vbuf_[j] = e.idx[j] >= 0 ? Vn_[e.idx[j]] : 0.0;
    return vbuf_;
}

bool Network::checkNonlinear(double t, Valves valves) {
    bool changed = false;
    for (std::size_t n : nonlinear_) {
        const Element& e = *comp_list_[n];
        if (e.valve() && (valves == NoValves || (valves == TurnOnOnly && e.valveOn()))) continue;
        changed = comp_list_[n]->check(t, terminalVoltages(*comp_list_[n])) || changed;
    }
    return changed;
}

void Network::resolve(double t) {
    refactor();
    calcBrnHistory(t); calcinjection(); Idivide(); calIUnew(); calcnewV(); reconstruct_V();
    ++resolves_;
}

void Network::calcNewbranchI(double t) {
    // Current delivered by the known sources into each known node: I = G_KU v_U + G_KK v_K - I_K
    Vector Iknown;
    if (src_count_ > 0) {
        const Vector a = G_KU_ * V_U_, b = G_KK_ * V_K_;
        Iknown.resize(src_count_);
        for (int k = 0; k < src_count_; ++k) Iknown[k] = a[k] + b[k] - I_K_[k];
    }
    for (std::size_t n = 0; n < comp_list_.size(); ++n) {
        auto& e = comp_list_[n];
        if (!knownIdx_[n].empty()) {
            Vector ik(knownIdx_[n].size());
            for (std::size_t j = 0; j < ik.size(); ++j) ik[j] = Iknown[knownIdx_[n][j]];
            e->setKnownCurrents(ik);
        }
        e->update(t, terminalVoltages(*e));
    }
}

// ---------------------------------------------------------------- chatter removal / interpolation

Vector Network::save() const {
    Vector out;
    out.reserve(stateSize_ ? stateSize_ : 4 * Vn_.size());
    out.insert(out.end(), Vn_.begin(), Vn_.end());
    for (const auto& e : comp_list_) e->saveState(out);
    stateSize_ = out.size();
    return out;
}

void Network::load(const Vector& s) {
    std::copy(s.begin(), s.begin() + numNodes_, Vn_.begin());
    std::size_t pos = numNodes_;
    for (auto& e : comp_list_) pos = e->loadState(s, pos);
}

void Network::solveStep(double t, Valves valves, int maxIterations) {
    Vdivide(); Vkupdate(t); calcBrnHistory(t); calcinjection(); Idivide(); calIUnew(); calcnewV(); reconstruct_V();
    int it = 0;
    while (checkNonlinear(t, valves)) {
        if (++it > maxIterations) { ++unconverged; break; }
        resolve(t);
    }
    calcNewbranchI(t);
}

double Network::detect(double a, double dt, bool interp, std::vector<Element*>& hits) const {
    double best = -1.0;
    hits.clear();
    for (Element* e : events_) {
        double al = e->eventAlpha(a, dt);
        if (al < 0.0) continue;
        if (!interp) al = 1.0;
        if (best < 0.0 || al < best - 1e-9) {
            best = al;
            hits.assign(1, e);
        } else if (al <= best + 1e-9) {
            hits.push_back(e);
        }
    }
    return best;
}

double Network::nodeVoltage(int node) const { return node == 0 ? 0.0 : Vn_[idx_.at(node)]; }

void Network::recordv(double t) {
    time_.push_back(t);
    Vector row(userNodes_);
    for (int n = 1; n <= userNodes_; ++n) row[n - 1] = nodeVoltage(n);
    trace_.push_back(std::move(row));
    Vector cur;
    cur.reserve(currentNames_.size());
    for (const auto& e : comp_list_) e->currents(cur);
    currentTrace_.push_back(std::move(cur));
}

void Network::writeResults(const std::string& path) const {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write " + path);
    f << "time";
    for (int n = 1; n <= userNodes_; ++n) f << ",V" << n;
    f << "\n" << std::setprecision(17);   // full double precision, so results compare exactly
    for (std::size_t s = 0; s < time_.size(); ++s) {
        f << time_[s];
        for (double v : trace_[s]) f << "," << v;
        f << "\n";
    }
}

void Network::writeCurrents(const std::string& path) const {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write " + path);
    f << "time";
    for (const auto& nm : currentNames_) f << "," << nm;
    f << "\n" << std::setprecision(17);
    for (std::size_t s = 0; s < time_.size(); ++s) {
        f << time_[s];
        for (double v : currentTrace_[s]) f << "," << v;
        f << "\n";
    }
}

void Network::writeEvents(const std::string& path) const {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write " + path);
    f << "type,t_start,t_end\n" << std::setprecision(12);
    for (const auto& e : comp_list_) {
        if (auto* fl = dynamic_cast<const Fault*>(e.get())) {
            f << "fault," << fl->t_apply << ",";
            if (fl->t_clear) f << *fl->t_clear;
            f << "\n";
        } else if (auto* sw = dynamic_cast<const Switch*>(e.get())) {
            for (double ts : sw->t_switch) f << "switch," << ts << ",\n";
        }
    }
}
