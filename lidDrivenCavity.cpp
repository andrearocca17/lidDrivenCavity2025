#include <iostream>
#include <fstream>
#include "Simple.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

int main() {
    // 1) Leggi config.json
    std::ifstream in("config.json");
    if (!in) {
        std::cerr << "ERROR: cannot open config.json\n";
        return 1;
    }
    json cfg; in >> cfg;

    // 2) Estrai parametri (default = valori originali)
    int    NX      = cfg.value("NX",      50);
    int    NY      = cfg.value("NY",      50);
    double L       = cfg.value("L",       1.0);
    double tol     = cfg.value("tol",     1e-6);
    int    maxIter = cfg.value("maxIter", 2000);
    int    outInt  = cfg.value("outInt",  200);

    std::cout << "[main] Config: "
              << "NX="      << NX
              << " NY="     << NY
              << " L="      << L
              << " tol="    << tol
              << " maxIter="<< maxIter
              << " outInt=" << outInt
              << "\n";

    // 3) Costruisci e lancia il solver **senza toccare il physics**
    Simple solver(NX, NY, L);
    solver.setTol(tol);
    solver.setMaxIter(maxIter);
    solver.setOutInt(outInt);
    solver.assembleSolveMomentum();

    return 0;
}
