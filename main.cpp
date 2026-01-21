#include <vector>
#include <fstream>
#include <array>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <ceres/ceres.h>

#include <nlohmann/json.hpp>

#include "model.h"

using json = nlohmann::json;

struct ResidualFunctor final : public ceres::CostFunction {
    ResidualFunctor(
        double t0_exp,              // Time of first experimental point
        double dt_exp,              // Experimental time increment
        std::vector<double> X_exp   // Experimental concentrations
    )
        : t0_exp_(t0_exp)
        , dt_exp_(dt_exp)
        , X_exp_(std::move(X_exp))
    {
        const int M = static_cast<int>(X_exp_.size());

        set_num_residuals(M);
        mutable_parameter_block_sizes()->push_back(5);
    }

    bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const override {

        const double* theta = parameters[0];
        const int M = static_cast<int>(X_exp_.size());

        Model::FittedParameters fitted_parameters {
            theta[0], theta[1], theta[2], theta[3], theta[4]
        };

        Model::FixedParameters fixed_parameters {
            .Di = 40.0,
            .R = 1.56 / 2.0,
            .L = 2.0,
            .F = 2.493333333333333 * 760.0 / 1.965,
            .X_in = 29424160260.695,
            .t_ads_start =  266.3915963445649,
            .t_ads_end = 543.5885793335225,
            .k_ads_smooth = 2.645834006658198,
            .dt = dt_exp_
        };

        Model model(fixed_parameters, fitted_parameters);

        // Skip until the first experimental point
        while (model.get_t() <= t0_exp_)
            model.do_step();

        std::vector<double> X_model(M);
        std::vector<std::array<double,5>> dX(M);

        for (long i = 0; i < M; i ++) {
            X_model[i] = model.get_Y()[0];

            for (long j = 0; j < 5; j ++)
                dX[i][j] = model.get_Y()[Model::sbase(j)];

            model.do_step();
        }

        // Residuals: r_i = X_model(t_i) - X_obs_i
        for (int i = 0; i < M; ++i) {
            residuals[i] = X_model[i] - X_exp_[i];
        }

        if (jacobians && jacobians[0]) {
            double* J = jacobians[0];
            // J[i*5 + j] = dr_i/dtheta_j = dX_i/dtheta_j
            for (int i = 0; i < M; ++i) {
                for (int j = 0; j < 5; ++j) {
                    J[i*5 + j] = dX[i][j];
                }
            }
        }

        return true;
    }

private:
    const double t0_exp_;
    const double dt_exp_;
    const std::vector<double> X_exp_;
};

int main() {

    std::ifstream ifs("../experimental_data.json");

    if (!ifs.good()) {
        fmt::println(stderr, "Error reading experimental data");
        exit(EXIT_FAILURE);
    }

    json data = json::parse(ifs);


    auto t_exp = data.at("t_exp").get<std::vector<double>>();
    auto X_exp = data.at("X_exp").get<std::vector<double>>();

    double t0_exp = t_exp[0];
    double dt_exp = 0.0;
    for (long i = 1; i < t_exp.size(); i ++) {
        dt_exp += t_exp[i] - t_exp[i-1];
    }
    dt_exp /= static_cast<double>(t_exp.size() - 1);

    auto * cost = new ResidualFunctor(t0_exp, dt_exp, X_exp);

    // Initial guesses
    double theta[] = {
            2.628413090171913e-12,  // k_ads
            0.02358734401037852,    // k_des
            1.3703757029773324e-16, // k_rxn
            58264568110461.61,      // S0
            96507325673713.16       // Y0
    };

//    Model::FittedParameters fitted_parameters_0 {
//            theta[0], theta[1], theta[2], theta[3], theta[4]
//    };
//
//    Model::FixedParameters fixed_parameters {
//            .Di = 40.0,
//            .R = 1.56 / 2.0,
//            .L = 2.0,
//            .F = 2.493333333333333 * 760.0 / 1.965,
//            .X_in = 29424160260.695,
//            .t_ads_start =  266.3915963445649,
//            .t_ads_end = 543.5885793335225,
//            .k_ads_smooth = 2.645834006658198,
//            .dt = dt_exp
//    };

//    Model model(fixed_parameters, fitted_parameters_0);

    ceres::Problem problem;
    problem.AddResidualBlock(cost, nullptr, theta);

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    // If your EvaluateAll is not thread-safe (CVODE usually isn't), force 1 thread:
    options.num_threads = 1;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";

    for (long i = 0; i < 5; i ++) {
        fmt::println("{}", theta[i]);
    }

    // Model::FixedParameters fixed_params{
    //     .Di = 40.0,
    //     .R = 1.56 / 2.0,
    //     .L = 2.0,
    //     .F = 2.493333333333333 * 760.0 / 1.965,
    //     .X_in = 29424160260.695,
    //     .t_ads_start =  266.3915963445649,
    //     .t_ads_end = 543.5885793335225,
    //     .k_ads_smooth = 2.645834006658198,
    //     .dt = 1.0
    // };
    //
    // Model::FittedParameters fitted_params{
    //     .k_ads  = 3.827187367718609e-12,
    //     .k_des  = 0.04253902879091293,
    //     .k_rxn  = 2.0493834170301142e-16,
    //     .S_tot  = 34482079719017.1,
    //     .P_tot  = 88448995701717.02
    // };
    //
    // Model model(fixed_params, fitted_params);
    //
    // const long N_t = 800;
    // std::vector<double> ts(N_t), Xs(N_t);
    // std::vector<std::array<double, 5>> dXs(N_t, std::array<double, 5>{});
    //
    // ts[0] = 0.0;
    // Xs[0] = fixed_params.X_in;
    // for (long i = 1; i < N_t; i ++) {
    //
    //     model.do_step();
    //
    //     ts[i] = model.get_t();
    //     Xs[i] = model.get_Y()[0];
    //     for (int j = 0; j < 5; j ++)
    //         dXs[i][j] = model.get_Y()[Model::sbase(j)];
    // }
    //
    // std::ofstream ofs("out_data.csv");
    //
    // fmt::println(ofs, "t, X, dX/dk_ads, dX/dk_des, dX/dk_rxn, dX/dS0, dX/dY0");
    //
    // double multipliers[] = {fitted_params.k_ads, fitted_params.k_des, fitted_params.k_rxn, fitted_params.S_tot, fitted_params.P_tot};
    //
    // for (long i = 0; i < N_t; i ++) {
    //     fmt::print(ofs, "{}, {}", ts[i], Xs[i]);
    //     for (int j = 0; j < 5; j ++) {
    //         fmt::print(ofs, ", {}", dXs[i][j]*multipliers[j]);
    //     }
    //     fmt::println(ofs, "");
    // }

    return 0;
}