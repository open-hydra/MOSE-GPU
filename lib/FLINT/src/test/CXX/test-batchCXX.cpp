#include "cantera/zerodim.h"
#include "cantera/numerics/Integrator.h"
#include "utils.h"
#include <vector>
#include <fstream>
#include <utility>
#include "cantera/base/AnyMap.h"
#include <string>
#include "cantera/base/global.h"


using namespace Cantera;
using std::cout;
using std::endl;
#include <chrono>

std::vector<std::pair<double, double>> simulateReactorTemperature(ReactorNet& sim, ThermoPhase& gas, int nsteps, double dt)
{
    std::vector<std::pair<double, double>> timeTemp;

    for (int i = 0; i <= nsteps; i++) {
        double tm = i * dt;
        sim.advance(tm);

        // AnyMap stats = sim.solverStats();
        // std::cout << "\n=== CVODE Solver Stats ===\n";
        // std::cout << "Steps taken:                 " << stats["steps"].asInt() << "\n";
        // std::cout << "RHS evaluations:             " << stats["rhs_evals"].asInt() << "\n";
        // std::cout << "Nonlinear iterations:        " << stats["nonlinear_iters"].asInt() << "\n";
        // std::cout << "Nonlinear convergence fails: " << stats["nonlinear_conv_fails"].asInt() << "\n";
        // std::cout << "Error test failures:         " << stats["err_test_fails"].asInt() << "\n";
        // std::cout << "Last method order used:      " << stats["last_order"].asInt() << "\n";
        // std::cout << "Stability limit reductions:  " << stats["stab_order_reductions"].asInt() << "\n";

        // std::cout << "Jacobian evaluations:        " << stats["jac_evals"].asInt() << "\n";
        // std::cout << "Linear solver setups:        " << stats["lin_solve_setups"].asInt() << "\n";
        // std::cout << "Linear RHS evaluations:      " << stats["lin_rhs_evals"].asInt() << "\n";
        // std::cout << "Linear iterations:           " << stats["lin_iters"].asInt() << "\n";
        // std::cout << "Linear convergence failures: " << stats["lin_conv_fails"].asInt() << "\n";

        // std::cout << "Preconditioner evaluations:  " << stats["prec_evals"].asInt() << "\n";
        // std::cout << "Preconditioner solves:       " << stats["prec_solves"].asInt() << "\n";
        // std::cout << "Jt-vector setup evals:       " << stats["jt_vec_setup_evals"].asInt() << "\n";
        // std::cout << "Jt-vector product evals:     " << stats["jt_vec_prod_evals"].asInt() << "\n";

        // if (stats.hasKey("step_solve_fails")) {
        //     std::cout << "Step solve failures:         " << stats["step_solve_fails"].asInt() << "\n";
        // }
        
        double T = gas.temperature();
        timeTemp.emplace_back(tm, T);
    }

    return timeTemp;
}

int main()
{
    try {

        std::vector<std::pair<std::string, double>> summaryData;

        std::chrono::high_resolution_clock::time_point t0;
        std::chrono::high_resolution_clock::time_point t1;
        std::chrono::duration<double> elapsed = t1 - t0;

        int nsteps = 0;
        double dt = 0;
        int sim_type = 0;

        std::cout << "What kind of simulation do you want to run?\n";
        std::cout << "1) verification\n";
        std::cout << "2) performance\n";
        std::cin >> sim_type;

        if (sim_type == 1) {
            nsteps = 1000;
        } else if (sim_type == 2) {
            nsteps = 1;
        } else {
            std::cerr << "Choose 1 or 2!\n";
            exit(EXIT_FAILURE);
        }

        std::cout << "Running with nstep = " << nsteps << "\n";

        //-------------------------------------------------------------------------------------------------
        // WD
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol1 = newSolution("WD/INPUT/WD.yaml");
        auto& gas = *sol1->thermo();
        gas.setState_TPY(1000.0, 100000.0, "CH4:0.20, O2:0.8");

        // Reactor setup
        IdealGasReactor reactor(sol1);
        ReactorNet sim;
        sim.addReactor(reactor);
        //sim.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 8.e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp1 = simulateReactorTemperature(sim, gas, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "WD Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("WD", elapsed.count());

        // Write results for first simulation
        std::ofstream outFile("WD/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp1) {
            outFile << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Troyes
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol2 = newSolution("Troyes/INPUT/troyes.yaml");
        auto& gas2 = *sol2->thermo();
        gas2.setState_TPY(1000.0, 100000.0, "H2:0.00534, CL2:0.18796, N2:0.80670");

        // Reactor setup
        IdealGasReactor reactor2(sol2);
        ReactorNet sim2;
        sim2.addReactor(reactor2);
        sim2.setTolerances(1e-12, 1e-15);

        // Integration settings
        dt = 5.0e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp2 = simulateReactorTemperature(sim2, gas2, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "Troyes Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("Troyes", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile2("Troyes/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp2) {
            outFile2 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Ecker
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol3 = newSolution("Ecker/INPUT/ecker.yaml");
        auto& gas3 = *sol3->thermo();
        gas3.setState_TPY(1000.0, 100000.0, "H2:0.00534534, CL2:0.18798856, N2:0.8066661");

        // Reactor setup
        IdealGasReactor reactor3(sol3);
        ReactorNet sim3;
        sim3.addReactor(reactor3);
        sim3.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 2.0e-2 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp3 = simulateReactorTemperature(sim3, gas3, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "Ecker Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("Ecker", elapsed.count());


        // Write results for second simulation
        std::ofstream outFile3("Ecker/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp3) {
            outFile3 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Cross
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol4 = newSolution("Cross/INPUT/cross.yaml");
        auto& gas4 = *sol4->thermo();
        gas4.setState_TPY(1010.0, 100000.0, "H2:0.00534534, CL2:0.18798856, N2:0.8066661");

        // Reactor setup
        IdealGasReactor reactor4(sol4);
        ReactorNet sim4;
        sim4.addReactor(reactor4);
        sim4.setTolerances(1e-7, 1e-7);
        sim4.setMaxSteps(100000);

        // Integration settings
        dt = 1.0e-2 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp4 = simulateReactorTemperature(sim4, gas4, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "Cross Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("Cross", elapsed.count());


        // Write results for second simulation
        std::ofstream outFile4("Cross/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp4) {
            outFile4 << time << "\t" << temp << "\n";
        }


        //-------------------------------------------------------------------------------------------------
        // Smooke
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol5 = newSolution("Smooke/INPUT/smooke.yaml");
        auto& gas5 = *sol5->thermo();
        gas5.setState_TPY(1300.0, 100000.0, "CH4:0.0552, O2:0.2201, N2: 0.7247");

        // Reactor setup
        IdealGasReactor reactor5(sol5);
        ReactorNet sim5;
        sim5.addReactor(reactor5);
        sim5.setTolerances(1e-7, 1e-7);
        sim5.setMaxSteps(1000000);

        // Integration settings
        dt = 0.2 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp5 = simulateReactorTemperature(sim5, gas5, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "Smooke Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("Smooke", elapsed.count());


        // Write results for second simulation
        std::ofstream outFile5("Smooke/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp5) {
            outFile5 << time << "\t" << temp << "\n";
        }


        //-------------------------------------------------------------------------------------------------
        // CORIA-CNRS
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol6 = newSolution("CORIA/INPUT/coria.yaml");
        auto& gas6 = *sol6->thermo();
        gas6.setState_TPY(1300.0, 100000.0, "CH4:0.2, O2:0.8");

        // Reactor setup
        IdealGasReactor reactor6(sol6);
        ReactorNet sim6;
        sim6.addReactor(reactor6);
        sim6.setTolerances(1e-7, 1e-7);
        sim6.setMaxSteps(1000000);

        // Integration settings
        dt = 0.005 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp6 = simulateReactorTemperature(sim6, gas6, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "CORIA Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("CORIA", elapsed.count());


        // Write results for second simulation
        std::ofstream outFile6("CORIA/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp6) {
            outFile6 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // TSR-CDF-13
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol7 = newSolution("TSR-CDF-13/INPUT/TSR-CDF-13.yaml");
        auto& gas7 = *sol7->thermo();
        gas7.setState_TPY(1300.0, 500000, "CH4:1, O2:4");

        // Reactor setup
        IdealGasReactor reactor7(sol7);
        ReactorNet sim7;
        sim7.addReactor(reactor7);
        sim7.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 5.0e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp7 = simulateReactorTemperature(sim7, gas7, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "TSR-CDF-13 Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("TSR-CDF-13", elapsed.count());


        // Write results for second simulation
        std::ofstream outFile7("TSR-CDF-13/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp7) {
            outFile7 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Pelucchi
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol8 = newSolution("Pelucchi/INPUT/pelucchi.yaml");
        auto& gas8 = *sol8->thermo();
        gas8.setState_TPY(1250.0, 100000, "CO:0.00859, O2:0.00606, H2O:0.00365, HCL:0.00025, N2:0.98044");
        //gas7.setState_TPX(1250.0, OneAtm, "CO:0.0086, O2:0.0053, H2O:0.0057, HCL:0.00019, N2:0.98021");

        // Reactor setup
        IdealGasReactor reactor8(sol8);
        ReactorNet sim8;
        sim8.addReactor(reactor8);
        sim8.setTolerances(1e-12, 1e-15);

        // Integration settings
        dt = 5.0e-2 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp8 = simulateReactorTemperature(sim8, gas8, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "Pelucchi Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("Pelucchi", elapsed.count());


        // Write results for second simulation
        std::ofstream outFile8("Pelucchi/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp8) {
            outFile8 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // ZK
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol9 = newSolution("ZK/INPUT/ZK.yaml");
        auto& gas9 = *sol9->thermo();
        gas9.setState_TPY(1300.0, 500000, "CH4:1, O2:4");

        // Reactor setup
        IdealGasReactor reactor9(sol9);
        ReactorNet sim9;
        sim9.addReactor(reactor9);
        sim9.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 5.0e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp9 = simulateReactorTemperature(sim9, gas9, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "ZK Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("ZK", elapsed.count());


        // Write results for second simulation
        std::ofstream outFile9("ZK/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp9) {
            outFile9 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // TSR-GP-24
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol10 = newSolution("TSR-GP-24/INPUT/TSR-GP-24.yaml");
        auto& gas10 = *sol10->thermo();
        gas10.setState_TPY(1300.0, 500000, "CH4:1, O2:4");

        // Reactor setup
        IdealGasReactor reactor10(sol10);
        ReactorNet sim10;
        sim10.addReactor(reactor10);
        sim10.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 5.0e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp10 = simulateReactorTemperature(sim10, gas10, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "TSR-GP-24 Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("TSR-GP-24", elapsed.count());


        // Write results for second simulation
        std::ofstream outFile10("TSR-GP-24/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp10) {
            outFile10 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // TSR-Rich-31
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol11 = newSolution("TSR-Rich-31/INPUT/TSR-Rich-31.yaml");
        auto& gas11 = *sol11->thermo();
        gas11.setState_TPY(1300.0, 500000, "CH4:1, O2:4");

        // Reactor setup
        IdealGasReactor reactor11(sol11);
        ReactorNet sim11;
        sim11.addReactor(reactor11);
        sim11.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 5.0e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp11 = simulateReactorTemperature(sim11, gas11, nsteps, dt);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "TSR-Rich-31 Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;
        summaryData.emplace_back("TSR-Rich-31", elapsed.count());


        // Write results for second simulation
        std::ofstream outFile11("TSR-Rich-31/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp11) {
            outFile11 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Summary output
        //-------------------------------------------------------------------------------------------------


        std::ofstream summaryFile("comp-batch-cantera.dat");
        for (const auto& [nspec, time] : summaryData) {
            summaryFile << nspec << "\t" << std::scientific << time << "\n";
        }

        appdelete();
        return 0;

    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        appdelete();
        return -1;
    }
}
