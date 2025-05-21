#include "Simple.h"
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

int main() {
    // 1) Read config.json
    std::ifstream in("config.json");
    if (!in) {
        std::cerr << "ERROR: cannot open config.json\n";
        return 1;
    }
    json cfg;
    in >> cfg;

    // 2) Extract mesh & domain (defaults to your old values)
    int N_    = cfg.value("NX",  5);
    int M_    = cfg.value("NY",  5);
    double L_ = cfg.value("L",   1.0);

    std::cout << "Using NX=" << N_ 
              << " NY=" << M_ 
              << " L="  << L_ << "\n";

    // 3) Call your existing solver as before
    Simple* sim = new Simple(N_, M_, L_);
    sim->assembleSolveMomentum();
    delete sim;
    return 0;
}
