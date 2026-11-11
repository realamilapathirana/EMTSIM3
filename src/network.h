#pragma once

// The network solver: Dommel's nodal method.
//
// Every time step the network is a resistive circuit  G v = i_h.  Each element
// adds its conductance block (stamp) and history currents. Nodes are split into
// unknown (U) and known (K, fixed by ideal sources and machine EMFs) nodes:
//
//     G_UU v_U = i_U - G_UK v_K
//
// G is stored as a sparse matrix. G_UU is symmetric positive definite, so it is
// factorised with a sparse LDL^T (minimum-degree ordering, no pivoting), at the
// start and again only when a switch operates or a non-linear element changes
// state. Non-linear elements (diodes, arresters, saturation) are handled by
// re-solving the step until their states agree with the solution.

#include <memory>
#include <string>
#include <vector>

#include "elements.h"
#include "linalg.h"

class Network {
public:
    // userNodes: nodes 1..userNodes are recorded; higher numbers are internal nodes of
    // composite elements (-1 = record all)
    explicit Network(std::vector<std::unique_ptr<Element>> elements, int userNodes = -1);

    bool updateSwitches(double t);  // scheduled switching; true if any conductance changed

    void formG();            // assemble the conductance matrix
    void Gdivide();          // split G into G_UU, G_UK, G_KU, G_KK
    void LU();               // factorise G_UU (sparse LDL^T)
    void refactor() { formG(); Gdivide(); LU(); }
    void Vdivide();          // split node voltages into V_U, V_K
    void Vkupdate(double t); // known node voltages V_K (ideal sources, machine EMFs)
    void calcBrnHistory(double t);   // history currents of every element
    void calcinjection();    // nodal current injections I_H
    void Idivide();          // split I_H into I_U, I_K
    void calIUnew();         // I_U_new = I_U - G_UK * V_K
    void calcnewV();         // solve G_UU V_U = I_U_new
    void reconstruct_V();    // Vn = [V_U ; V_K]
    enum Valves { NoValves, TurnOnOnly, AllValves };     // which diodes/thyristors check() may switch
    bool checkNonlinear(double t, Valves valves = AllValves);   // state changes of non-linear elements; true -> re-solve
    void resolve(double t);  // refactorise and solve the step again (after checkNonlinear)
    void calcNewbranchI(double t);   // commit: element currents and states for the next step
    void recordv(double t);  // record node voltages and branch currents

    // chatter removal / switching interpolation
    Vector save() const;                 // continuous state: node voltages and every element's state
    void load(const Vector& s);
    void solveStep(double t, Valves valves, int maxIterations);   // one step with the present switch states
    double detect(double a, double dt, bool interp, std::vector<Element*>& hits) const;   // < 0: none
    void runWithEvents(int maxIterations);   // the time loop with both (config::chatter / config::interp)
    int unconverged = 0;                 // steps where non-linear elements did not settle

    void writeResults(const std::string& path) const;   // time, V1..Vn as CSV
    void writeCurrents(const std::string& path) const;  // time, branch currents as CSV
    void writeEvents(const std::string& path) const;    // fault/switch times as CSV

    int numNodes() const { return userNodes_; }        // recorded (user) nodes
    int totalNodes() const { return numNodes_; }       // including internal nodes
    double nodeVoltage(int node) const;                // node number, 0 = ground
    const std::vector<double>& time() const { return time_; }
    const std::vector<Vector>& trace() const { return trace_; }
    const std::vector<std::string>& currentNames() const { return currentNames_; }
    const std::vector<Vector>& currentTrace() const { return currentTrace_; }
    const std::vector<std::unique_ptr<Element>>& elements() const { return comp_list_; }
    std::size_t nnzG_UU() const { return G_UU_.nnz(); }                 // non-zeros in G_UU
    std::size_t nnzFactor() const { return G_UU_LU_.nnzL(); }           // non-zeros in L (below the diagonal)
    int factorisations() const { return factorisations_; }
    int resolves() const { return resolves_; }

private:
    const Vector& terminalVoltages(const Element& e);   // into a reused buffer

    std::vector<std::unique_ptr<Element>> comp_list_;
    int numNodes_ = 0, userNodes_ = 0, unknownNodes_ = 0, src_count_ = 0;
    std::vector<int> order_;  // internal index -> node
    std::vector<int> idx_;    // node -> internal index (ground -> -1)

    SparseMatrix G_, G_UU_, G_UK_, G_KU_, G_KK_;
    SparseLDL G_UU_LU_;
    int factorisations_ = 0, resolves_ = 0;
    Vector Vn_, V_U_, V_K_, I_H_, I_U_, I_K_, I_U_new_;
    std::vector<Vector> hist_;                  // history currents of each element
    std::vector<std::vector<int>> knownIdx_;    // per element: indices into V_K of its known terminals
    std::vector<std::size_t> nonlinear_;        // elements whose state can change with the solution
    std::vector<Element*> events_;              // switches, faults, diodes, thyristors
    mutable std::size_t stateSize_ = 0;         // length of save()'s vector, once known
    Vector vbuf_;                               // terminal voltages of one element

    std::vector<double> time_;
    std::vector<Vector> trace_;         // one row per time step: voltages of nodes 1..userNodes
    std::vector<std::string> currentNames_;
    std::vector<Vector> currentTrace_;  // one row per time step: all branch currents
};
