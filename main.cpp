#include <vector>
#include <fstream>
#include <array>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <ceres/ceres.h>

#include <nlohmann/json.hpp>

#include "model.h"

using json = nlohmann::json;
using ordered_json = nlohmann::ordered_json;

struct ResidualFunctor final : public ceres::CostFunction {
    ResidualFunctor(
        double t0_exp,                                      // Time of first experimental point
        Model::FixedParameters const & fixed_parameters,    // Experimental parameters that are fixed
        std::vector<double> X_exp                           // Experimental concentrations
    )
        : t0_exp_(t0_exp)
        , fixed_parameters_(fixed_parameters)
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

        Model model(fixed_parameters_, fitted_parameters);

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
    const Model::FixedParameters fixed_parameters_;
    const std::vector<double> X_exp_;
};

struct InputData {
    double t0_exp;
    Model::FixedParameters fixed_parameters;
    Model::FittedParameters initial_guess;
    std::vector<double> t_exp, X_exp;
};

InputData load_input_data(std::filesystem::path const & input_file_path) {

    std::ifstream ifs(input_file_path);

    if (!ifs.good()) {
        fmt::println(stderr, "Error reading experimental data from {}", input_file_path.string());
        exit(EXIT_FAILURE);
    }

    json data = json::parse(ifs);

    auto t_exp = data.at("experimental_data").at("t_exp").get<std::vector<double>>();

    double t0_exp = t_exp[0];
    double dt_exp = 0.0;
    for (long i = 1; i < t_exp.size(); i ++) {
        dt_exp += t_exp[i] - t_exp[i-1];
    }
    dt_exp /= static_cast<double>(t_exp.size() - 1);

    InputData input_data{
        .t0_exp = t0_exp,
        .fixed_parameters = {
            .Di = 40.0,
            .R = data.at("R").get<double>(),
            .L = data.at("L").get<double>(),
            .F = data.at("F").get<double>() * 760.0 / data.at("pressure").get<double>(),
            .X_in = data.at("X_in").get<double>(),
            .t_ads_start = data.at("t_ads_start").get<double>(),
            .t_ads_end = data.at("t_ads_end").get<double>(),
            .k_ads_smooth = 2.0,
            .dt = dt_exp
        },
        .initial_guess = {
            .k_ads = data.at("initial_guess").at("k_ads").get<double>(),
            .k_des = data.at("initial_guess").at("k_des").get<double>(),
            .k_rxn = data.at("initial_guess").at("k_rxn").get<double>(),
            .S_tot = data.at("initial_guess").at("S_tot").get<double>(),
            .P_tot = data.at("initial_guess").at("Y_tot").get<double>()
        },
        .t_exp = data.at("experimental_data").at("t_exp").get<std::vector<double>>(),
        .X_exp = data.at("experimental_data").at("X_exp").get<std::vector<double>>()
    };

    return input_data;
}

std::pair<std::vector<double>, std::vector<std::array<double, 5>>> solve_model(
    Model::FixedParameters const & fixed_parameters,
    Model::FittedParameters const & fitted_parameters,
    const long N_t) {

    Model model(fixed_parameters, fitted_parameters);

    std::vector<double> ts(N_t), Xs(N_t);
    std::vector<std::array<double, 5>> dXs(N_t, std::array<double, 5>{});

    ts[0] = 0.0;
    Xs[0] = fixed_parameters.X_in;
    for (long i = 1; i < N_t; i ++) {

        model.do_step();

        ts[i] = model.get_t();
        Xs[i] = model.get_Y()[0];
        for (int j = 0; j < 5; j ++)
            dXs[i][j] = model.get_Y()[Model::sbase(j)];
    }

    return std::make_pair(Xs, dXs);
}

int main(int argc, char ** argv) {

    if (argc < 2) {
        fmt::println(stderr, "Input file path must be provided as an argument");
        return EXIT_FAILURE;
    }

    std::filesystem::path input_file_path(argv[1]);
    std::filesystem::path output_file_path = input_file_path.parent_path() / "fitted.json";

    auto input_data = load_input_data(input_file_path);

    auto * cost = new ResidualFunctor(input_data.t0_exp, input_data.fixed_parameters, input_data.X_exp);

    // Initial guesses
    double theta[] = {
            input_data.initial_guess.k_ads,     // k_ads
            input_data.initial_guess.k_des,     // k_des
            input_data.initial_guess.k_rxn,     // k_rxn
            input_data.initial_guess.S_tot,     // S0
            input_data.initial_guess.P_tot      // Y0
    };

    Model::FittedParameters fitted_parameters_0 {
            theta[0], theta[1], theta[2], theta[3], theta[4]
    };

    auto [X0, dX0] = solve_model(input_data.fixed_parameters, fitted_parameters_0, input_data.t_exp.size());

    ceres::Problem problem;
    problem.AddResidualBlock(cost, nullptr, theta);

    // Add constraints
    problem.SetParameterLowerBound(theta, 0, 0.0);
    problem.SetParameterLowerBound(theta, 1, 0.0);
    problem.SetParameterLowerBound(theta, 2, 0.0);
    problem.SetParameterLowerBound(theta, 3, 0.0);
    problem.SetParameterLowerBound(theta, 4, 0.0);

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    // If your EvaluateAll is not thread-safe (CVODE usually isn't), force 1 thread:
    options.num_threads = 1;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";

    Model::FittedParameters fitted_parameters {
        theta[0], theta[1], theta[2], theta[3], theta[4]
};

    auto [X, dX] = solve_model(input_data.fixed_parameters, fitted_parameters, input_data.t_exp.size());

    ordered_json out_data;
    out_data["solution"]["k_ads"] = theta[0];
    out_data["solution"]["k_des"] = theta[1];
    out_data["solution"]["k_rxn"] = theta[2];
    out_data["solution"]["S_tot"] = theta[3];
    out_data["solution"]["Y_tot"] = theta[4];
    out_data["fitted_data"]["X0"] = X0;
    out_data["fitted_data"]["X"] = X;

    std::ofstream ofs(output_file_path);

    if (!ofs.good()) {
        fmt::println(stderr, "Error writing output to {}", output_file_path.string());
        exit(EXIT_FAILURE);
    }

    ofs << out_data.dump(4);

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