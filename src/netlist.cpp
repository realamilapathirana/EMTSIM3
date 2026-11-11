#include "netlist.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "config.h"

namespace {

constexpr double PI = 3.14159265358979323846;

// Number of node terminals of each element type
int terminalCount(const std::string& kind) {
    static const std::map<std::string, int> n = {
        {"R", 2}, {"L", 2}, {"C", 2}, {"V", 2}, {"SW", 2}, {"F", 2}, {"PI", 2}, {"TL", 2}, {"LSAT", 2},
        {"MOV", 2}, {"D", 2}, {"THY", 2}, {"XF", 4}, {"RL3", 6}, {"PI3", 6}, {"TL3", 6}, {"XF3", 6},
        {"V3", 3}, {"SM", 3}, {"IM", 3}};
    const auto it = n.find(kind);
    return it == n.end() ? -1 : it->second;
}

struct Statement {
    std::string kind, name, where;
    std::vector<int> nodes;
    std::vector<std::string> args;          // positional arguments after the nodes
    std::map<std::string, double> kw;       // key=value
    std::set<std::string> flags;            // bare words (czero)
};

std::string upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

double toNumber(const std::string& s, const std::string& what) {
    try {
        std::size_t pos = 0;
        const double v = std::stod(s, &pos);
        if (pos != s.size()) throw std::invalid_argument("");
        return v;
    } catch (const std::logic_error&) {
        throw std::runtime_error(what + " \"" + s + "\" is not a number");
    }
}

class Builder {
public:
    Builder(Circuit& c, int maxNode) : ckt_(c), next_(maxNode + 1) {}

    int internalNode() { return next_++; }

    template <class T, class... Args>
    T& add(const std::string& name, Args&&... args) {
        T& e = ckt_.add<T>(std::forward<Args>(args)...);
        std::string n = name;
        for (int k = 2; used_.count(n); ++k) n = name + "~" + std::to_string(k);
        used_.insert(n);
        e.name = n;
        return e;
    }

private:
    Circuit& ckt_;
    int next_;
    std::set<std::string> used_;
};

// Positional argument k (required) / optional with default; key=value overrides
struct Args {
    const Statement& s;
    double req(std::size_t k, const char* what) const {
        if (k >= s.args.size()) throw std::runtime_error(std::string("missing ") + what);
        return toNumber(s.args[k], what);
    }
    double opt(std::size_t k, const char* key, double def) const {
        const auto it = s.kw.find(key);
        if (it != s.kw.end()) return it->second;
        return k < s.args.size() ? toNumber(s.args[k], key) : def;
    }
    double kwOr(const char* key, double def) const {
        const auto it = s.kw.find(key);
        return it == s.kw.end() ? def : it->second;
    }
    bool has(std::size_t k) const { return k < s.args.size(); }
};

void build(const Statement& s, Builder& b) {
    const Args a{s};
    const auto& n = s.nodes;
    const std::string& name = s.name;
    const std::string& k = s.kind;

    if (k == "R" || k == "L" || k == "C") {
        b.add<RLC>(name, k[0], n[0], n[1], a.req(0, "value"));
    } else if (k == "V") {
        const double mag = a.req(0, "magnitude"), ang = a.req(1, "angle"), f = a.req(2, "frequency");
        const double R = a.opt(3, "r", 0.0), L = a.opt(4, "l", 0.0);
        if (L > 0.0) {   // ideal source behind series R + L
            if (n[0] != 0 && n[1] != 0) throw std::runtime_error("a source with inductance needs one terminal on ground");
            const int x = b.internalNode();
            Matrix Rm(1, 1), Lm(1, 1);
            Rm(0, 0) = R; Lm(0, 0) = L;
            if (n[0] == 0) {
                b.add<Source>(name, 0, x, mag, ang, f, 0.0);
                b.add<CoupledRL>(name + ".Z", std::vector<int>{x}, std::vector<int>{n[1]}, Rm, Lm);
            } else {
                b.add<Source>(name, x, 0, mag, ang, f, 0.0);
                b.add<CoupledRL>(name + ".Z", std::vector<int>{n[0]}, std::vector<int>{x}, Rm, Lm);
            }
        } else {
            b.add<Source>(name, n[0], n[1], mag, ang, f, R);
        }
    } else if (k == "SW") {
        if (s.args.empty() || (s.args[0] != "closed" && s.args[0] != "open"))
            throw std::runtime_error("SW state must be 'closed' or 'open'");
        std::vector<double> times;
        for (std::size_t i = 1; i < s.args.size(); ++i) times.push_back(toNumber(s.args[i], "switching time"));
        Switch& sw = b.add<Switch>(name, n[0], n[1], s.args[0] == "closed", times, a.kwOr("ron", 0.01), a.kwOr("roff", 1.0e6));
        sw.opens_at_current_zero = s.flags.count("czero") > 0 || s.kw.count("chop") > 0;
        sw.chop = a.kwOr("chop", 0.0);
    } else if (k == "F") {
        const double t0 = a.req(0, "t_apply");
        const std::optional<double> t1 = a.has(1) ? std::optional<double>(a.req(1, "t_clear")) : std::nullopt;
        b.add<Fault>(name, n[0], n[1], t0, t1, a.kwOr("ron", 0.01), a.kwOr("roff", 1.0e6));
    } else if (k == "RL3") {
        const double R1 = a.req(0, "R1"), L1 = a.req(1, "L1");
        auto& e = b.add<CoupledRL>(name, CoupledRL::fromSequence({n[0], n[1], n[2]}, {n[3], n[4], n[5]}, R1, L1,
                                                                 a.opt(2, "r0", R1), a.opt(3, "l0", L1)));
        e.type = "RL3";
    } else if (k == "PI") {
        const double R = a.req(0, "R"), L = a.req(1, "L"), Cc = a.req(2, "C");
        Matrix Rm(1, 1), Lm(1, 1);
        Rm(0, 0) = R; Lm(0, 0) = L;
        b.add<CoupledRL>(name, std::vector<int>{n[0]}, std::vector<int>{n[1]}, Rm, Lm).type = "PI";
        if (Cc > 0.0) {
            b.add<RLC>(name + ".Ck", 'C', n[0], 0, Cc / 2.0);
            b.add<RLC>(name + ".Cm", 'C', n[1], 0, Cc / 2.0);
        }
    } else if (k == "PI3") {
        const double R1 = a.req(0, "R1"), L1 = a.req(1, "L1"), C1 = a.req(2, "C1");
        const double R0 = a.opt(3, "r0", R1), L0 = a.opt(4, "l0", L1), C0 = a.opt(5, "c0", C1);
        b.add<CoupledRL>(name, CoupledRL::fromSequence({n[0], n[1], n[2]}, {n[3], n[4], n[5]}, R1, L1, R0, L0)).type = "PI3";
        if (C1 > 0.0) {
            b.add<CoupledC>(name + ".Ck", CoupledC::fromSequence({n[0], n[1], n[2]}, C1 / 2.0, C0 / 2.0));
            b.add<CoupledC>(name + ".Cm", CoupledC::fromSequence({n[3], n[4], n[5]}, C1 / 2.0, C0 / 2.0));
        }
    } else if (k == "V3") {
        const double mag = a.req(0, "magnitude"), ang = a.req(1, "angle"), f = a.req(2, "frequency");
        const double shift[3] = {0.0, -120.0, 120.0};
        const char* ph[3] = {".a", ".b", ".c"};
        const bool hasL = a.has(4) || s.kw.count("l1") > 0;     // an inductance is given
        const double R1 = a.opt(3, "r1", 0.0), L1 = a.opt(4, "l1", 0.0);
        if (hasL) {   // ideal EMFs behind a coupled (sequence) impedance
            std::vector<int> x;
            for (int p = 0; p < 3; ++p) {
                x.push_back(b.internalNode());
                b.add<Source>(name + ph[p], 0, x[p], mag, ang + shift[p], f, 0.0);
            }
            b.add<CoupledRL>(name + ".Z", CoupledRL::fromSequence(x, {n[0], n[1], n[2]}, R1, L1, a.opt(5, "r0", R1),
                                                                   a.opt(6, "l0", L1)));
        } else {
            for (int p = 0; p < 3; ++p) b.add<Source>(name + ph[p], 0, n[p], mag, ang + shift[p], f, R1);
        }
    } else if (k == "TL") {
        b.add<BergeronLine>(name, std::vector<int>{n[0]}, std::vector<int>{n[1]}, std::vector<double>{a.req(0, "Zc")},
                            std::vector<double>{a.req(1, "tau")}, std::vector<double>{a.opt(2, "r", 0.0)});
    } else if (k == "TL3") {
        const double Z1 = a.req(0, "Z1"), t1 = a.req(1, "tau1"), Z0 = a.req(2, "Z0"), t0 = a.req(3, "tau0");
        const double R1 = a.opt(4, "r1", 0.0), R0 = a.opt(5, "r0", R1);
        auto& e = b.add<BergeronLine>(name, std::vector<int>{n[0], n[1], n[2]}, std::vector<int>{n[3], n[4], n[5]},
                                      std::vector<double>{Z0, Z1, Z1}, std::vector<double>{t0, t1, t1},
                                      std::vector<double>{R0, R1, R1});
        e.type = "TL3";
    } else if (k == "XF" || k == "XF3") {
        const bool three = k == "XF3";
        std::string conn;
        std::size_t o = 0;
        if (three) {
            if (s.args.empty()) throw std::runtime_error("XF3 needs a connection, e.g. YNd1");
            conn = s.args[0];
            o = 1;
        }
        const double S = a.req(o + 0, "S"), V1 = a.req(o + 1, "V1"), V2 = a.req(o + 2, "V2"), f = a.req(o + 3, "f");
        const double x = a.req(o + 4, "x"), r = a.req(o + 5, "r");
        const double im = a.kwOr("im", 1.0), knee = a.kwOr("knee", 0.0), lsatPu = a.kwOr("lsat", 2.0 * x);
        const double w = 2.0 * PI * f;
        // windings: primary (from, to), secondary (from, to) per phase
        std::vector<std::array<int, 4>> wind;
        double Vw1 = V1, Vw2 = V2, Sw = S;
        std::vector<int> deltaNodes;
        if (!three) {
            wind.push_back({n[0], n[1], n[2], n[3]});
        } else {
            static const std::regex re("^(YN|Y|D)(yn|y|d)(0|1|11)$");
            std::smatch m;
            if (!std::regex_match(conn, m, re)) throw std::runtime_error("XF3 connection \"" + conn + "\": use e.g. YNyn0, YNd1, YNd11, Dyn1, Dyn11, Dd0");
            const std::string P = m[1], Q = m[2];
            const int clock = std::stoi(m[3]);
            const bool pd = P == "D", sd = Q == "d";
            if ((pd == sd && clock != 0) || (pd != sd && clock == 0))
                throw std::runtime_error("XF3 connection \"" + conn + "\": star-star and delta-delta need clock 0, star-delta 1 or 11");
            Vw1 = pd ? V1 : V1 / std::sqrt(3.0);
            Vw2 = sd ? V2 : V2 / std::sqrt(3.0);
            Sw = S / 3.0;
            const int Np = P == "Y" ? b.internalNode() : 0, Ns = Q == "y" ? b.internalNode() : 0;
            for (int p = 0; p < 3; ++p) {
                const int nxt = (p + 1) % 3, prv = (p + 2) % 3;
                // delta winding orientation sets the phase shift (clock 1: secondary lags 30 deg, 11: leads)
                const int pt = pd ? n[(sd || clock == 11) ? nxt : prv] : Np;
                const int st = sd ? n[3 + ((pd || clock == 1) ? nxt : prv)] : Ns;
                wind.push_back({n[p], pt, n[3 + p], st});
            }
            // a delta winding has no ground reference of its own: add a very high resistance to ground
            if (pd) deltaNodes.insert(deltaNodes.end(), {n[0], n[1], n[2]});
            if (sd) deltaNodes.insert(deltaNodes.end(), {n[3], n[4], n[5]});
        }
        const double Zb = Vw1 * Vw1 / Sw, nratio = Vw1 / Vw2;
        const char* ph[3] = {".a", ".b", ".c"};
        for (std::size_t p = 0; p < wind.size(); ++p) {
            const std::string nm = three ? name + ph[p] : name;
            auto& core = b.add<TransformerCore>(nm, wind[p][0], wind[p][1], wind[p][2], wind[p][3], nratio, r * Zb, x * Zb / w);
            core.type = k;
            if (im > 0.0) {
                const double Lm = Vw1 * Vw1 / (w * Sw * im / 100.0);
                if (knee > 0.0) {   // two-slope saturation curve: Lm up to the knee, then lsat
                    const double lamR = std::sqrt(2.0) * Vw1 / w, iK = knee * lamR / Lm, Ls = lsatPu * Zb / w;
                    b.add<SaturableInductor>(nm + ".m", wind[p][0], wind[p][1], std::vector<double>{iK, iK + lamR / Ls},
                                             std::vector<double>{knee * lamR, (knee + 1.0) * lamR});
                } else {
                    b.add<RLC>(nm + ".m", 'L', wind[p][0], wind[p][1], Lm);
                }
            }
        }
        for (int d : deltaNodes) b.add<RLC>(name + ".ref", 'R', d, 0, 1.0e8);
    } else if (k == "LSAT") {
        if (s.args.size() < 2 || s.args.size() % 2) throw std::runtime_error("LSAT needs pairs of current and flux");
        std::vector<double> i, l;
        for (std::size_t p = 0; p < s.args.size(); p += 2) {
            i.push_back(toNumber(s.args[p], "current"));
            l.push_back(toNumber(s.args[p + 1], "flux"));
        }
        b.add<SaturableInductor>(name, n[0], n[1], i, l);
    } else if (k == "MOV") {
        b.add<SurgeArrester>(name, n[0], n[1], a.req(0, "Vref"), a.opt(1, "iref", 1000.0), a.opt(2, "q", 26.0));
    } else if (k == "D") {
        b.add<Diode>(name, n[0], n[1], a.kwOr("ron", 1e-3), a.kwOr("roff", 1e6), a.kwOr("vf", 0.0));
    } else if (k == "THY") {
        b.add<Thyristor>(name, n[0], n[1], a.req(0, "alpha"), a.req(1, "frequency"), a.kwOr("ref", 0.0),
                         a.kwOr("width", 120.0), a.kwOr("ron", 1e-3), a.kwOr("roff", 1e6));
    } else if (k == "SM") {
        std::vector<int> x{b.internalNode(), b.internalNode(), b.internalNode()};
        b.add<SynchronousMachine>(name, std::vector<int>{n[0], n[1], n[2]}, x, a.req(0, "S"), a.req(1, "V"), a.req(2, "f"),
                                  a.req(3, "H"), a.req(4, "xd2"), a.req(5, "ra"), a.req(6, "E"), a.req(7, "delta0"),
                                  a.req(8, "Pm"), a.opt(9, "d", 0.0));
    } else if (k == "IM") {
        std::vector<int> x{b.internalNode(), b.internalNode(), b.internalNode()};
        b.add<InductionMachine>(name, std::vector<int>{n[0], n[1], n[2]}, x, a.req(0, "S"), a.req(1, "V"), a.req(2, "f"),
                                int(a.req(3, "poles")), a.req(4, "rs"), a.req(5, "xls"), a.req(6, "rr"), a.req(7, "xlr"),
                                a.req(8, "xm"), a.req(9, "H"), a.req(10, "Tload"), a.opt(11, "speed0", 0.0));
    } else {
        throw std::runtime_error("unknown element type '" + k + "'");
    }
}

}  // namespace

Circuit parseNetlist(const std::string& text, const std::string& where) {
    std::istringstream in(text);
    std::vector<Statement> stmts;
    std::string raw;
    int maxNode = 0;
    // time step and stop time: defaults unless the netlist has a ".tran dt tstop" line
    config::Dt = config::step = config::defaultDt;
    config::euler = false;
    config::stpTime = config::srtTime + config::defaultStop;
    config::chatter = config::defaultChatter;
    config::interp = config::defaultInterp;
    for (int lineno = 1; std::getline(in, raw); ++lineno) {
        Statement s;
        s.where = where + ":" + std::to_string(lineno);
        const std::size_t hash = raw.find('#');
        std::istringstream ss(raw.substr(0, hash));
        std::vector<std::string> tok;
        for (std::string t; ss >> t;) tok.push_back(t);
        if (tok.empty()) continue;
        if (upper(tok[0]) == ".TRAN") {
            if (tok.size() != 3) throw std::runtime_error(s.where + ": use .tran dt tstop");
            const double dt = toNumber(tok[1], "dt"), stop = toNumber(tok[2], "tstop");
            if (!(dt > 0.0) || !(stop > dt)) throw std::runtime_error(s.where + ": .tran needs 0 < dt < tstop");
            if (!stmts.empty()) throw std::runtime_error(s.where + ": put .tran before the elements");
            config::Dt = config::step = dt;
            config::stpTime = config::srtTime + stop;
            continue;
        }
        if (upper(tok[0]) == ".OPTIONS") {
            for (std::size_t k = 1; k < tok.size(); ++k) {
                const std::size_t eq = tok[k].find('=');
                const std::string key = lower(tok[k].substr(0, eq));
                const int on = eq == std::string::npos ? -1 : parseOnOff(tok[k].substr(eq + 1));
                if ((key != "chatter" && key != "interp") || on < 0)
                    throw std::runtime_error(s.where + ": use .options chatter=on|off interp=on|off");
                (key == "chatter" ? config::chatter : config::interp) = on;
            }
            continue;
        }
        if (hash != std::string::npos) {   // the first word of the comment names the element
            std::istringstream cs(raw.substr(hash + 1));
            cs >> s.name;
            while (!s.name.empty() && (s.name.back() == ',' || s.name.back() == ';' || s.name.back() == ':')) s.name.pop_back();
        }
        s.kind = upper(tok[0]);
        const int nt = terminalCount(s.kind);
        try {
            if (nt < 0) throw std::runtime_error("unknown element type '" + tok[0] + "'");
            if (int(tok.size()) < 1 + nt) throw std::runtime_error("expected " + std::to_string(nt) + " node numbers");
            for (int i = 0; i < nt; ++i) {
                const double v = toNumber(tok[1 + i], "node");
                if (v < 0 || v != std::floor(v)) throw std::runtime_error("node numbers must be whole numbers >= 0");
                s.nodes.push_back(int(v));
                maxNode = std::max(maxNode, int(v));
            }
            for (std::size_t i = 1 + nt; i < tok.size(); ++i) {
                const std::size_t eq = tok[i].find('=');
                if (eq != std::string::npos) {
                    const std::string key = lower(tok[i].substr(0, eq));
                    s.kw[key] = toNumber(tok[i].substr(eq + 1), key);
                } else if (tok[i] == "czero") {
                    s.flags.insert(tok[i]);
                } else {
                    s.args.push_back(tok[i]);
                }
            }
        } catch (const std::runtime_error& e) {
            throw std::runtime_error(s.where + ": " + e.what());
        }
        if (s.name.empty()) s.name = s.kind + "@" + std::to_string(lineno);
        stmts.push_back(std::move(s));
    }

    Circuit ckt;
    ckt.userNodes = maxNode;
    Builder b(ckt, maxNode);
    for (const auto& s : stmts) {
        try {
            build(s, b);
        } catch (const std::exception& e) {
            throw std::runtime_error(s.where + ": " + e.what());
        }
    }
    return ckt;
}

Circuit loadNetlist(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open netlist " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return parseNetlist(ss.str(), path);
}

int parseOnOff(const std::string& value) {
    const std::string v = lower(value);
    if (v == "on" || v == "1" || v == "true" || v == "yes") return 1;
    if (v == "off" || v == "0" || v == "false" || v == "no") return 0;
    return -1;
}
