#include <vector>
#include <fstream>
#include <array>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <ceres/ceres.h>

#include <nlohmann/json.hpp>

#include "model_chained.h"

// #define ENABLE_SCALING

#ifdef ENABLE_SCALING
#define S_EXP(__X__) exp(__X__)
#define S_LOG(__X__) log(__X__)
#else //ENABLE_SCALING
#define S_EXP(__X__) (__X__)
#define S_LOG(__X__) (__X__)
#endif //ENABLE_SCALING

using json = nlohmann::json;
using ordered_json = nlohmann::ordered_json;

bool almost_equal(double a, double b,
                  double rel_tol = 1e-2,
                  double abs_tol = 1e-2)
{
    return std::abs(a - b) <=
           std::max(rel_tol * std::max(std::abs(a), std::abs(b)),
                    abs_tol);
}

struct ResidualFunctor final : public ceres::CostFunction {
    ResidualFunctor(
        double t0_exp,                                      // Time of first experimental point
        Model::FixedParameters const & fixed_parameters,    // Experimental parameters that are fixed
        std::vector<double> X_exp,                          // Experimental concentrations
        long N_reactors
    )
        : t0_exp_(t0_exp)
        , fixed_parameters_(fixed_parameters)
        , X_exp_(std::move(X_exp))
        , N_reactors_(N_reactors)
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

        // fmt::println("{:.03e} {:.03e} {:.03e} {:.03e} {:.03e}", S_EXP(theta[0]), S_EXP(theta[1]), S_EXP(theta[2]), S_EXP(theta[3]), S_EXP(theta[4]));

        Model::FittedParameters fitted_parameters {
            S_EXP(theta[0]), S_EXP(theta[1]), S_EXP(theta[2]), S_EXP(theta[3]), S_EXP(theta[4])
        };

        Model model(fixed_parameters_, fitted_parameters, N_reactors_);

        // Skip until the first experimental point
        long alignment_steps = 100;
        while (!almost_equal(model.get_t(), t0_exp_) && alignment_steps > 0) {
            model.do_step();
            fmt::println("{} {}", t0_exp_, model.get_t());
            alignment_steps --;
        }

        if (alignment_steps <= 0) {
            fmt::println(stderr, "Failed to align model and experimental time grid");
            exit(EXIT_FAILURE);
        }

        std::vector<double> X_model(M);
        std::vector<std::array<double,5>> dX(M);

        for (long i = 0; i < M; i ++) {

            // fmt::println("{} {}", t0_exp_ + static_cast<double>(i) * 3.333333333333333, model.get_t());

            X_model[i] = model.get_Y()[model.xbase(N_reactors_-1, 0)];

            for (long j = 0; j < 5; j ++)
#ifndef ENABLE_SCALING
                dX[i][j] = model.get_Y()[model.sbase(N_reactors_-1, j)];
#else //ENABLE_SCALING
                dX[i][j] = exp(theta[j]) * model.get_Y()[model.sbase(N_reactors_-1, j)];
#endif //ENABLE_SCALING

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
    const long N_reactors_;
};

struct InputData {
    double t0_exp;
    Model::FixedParameters fixed_parameters;
    Model::FittedParameters initial_guess;
    std::vector<double> t_exp, X_exp;
    long N_reactors;
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

    fmt::println("dt_exp: {}", dt_exp);

    // Correct t_exp if not equally spaced (and correct t0_exp so it is a multiple of dt_exp)
    t0_exp = round(t0_exp / dt_exp) * dt_exp;
    for (long i = 0; i < t_exp.size(); i ++) {
        t_exp[i] = t0_exp + static_cast<double>(i) * dt_exp;
    }

    auto X_exp = data.at("experimental_data").at("X_exp").get<std::vector<double>>();

    // std::rotate(X_exp.rbegin(), X_exp.rbegin() + 1, X_exp.rend());
    // std::rotate(X_exp.begin(), X_exp.begin() + 1, X_exp.end());

    InputData input_data{
        .t0_exp = t0_exp,
        .fixed_parameters = {
            .Di = 40.0,
            .R = data.at("R").get<double>(),
            .L = data.at("L").get<double>(),
            .F = data.at("F").get<double>() * 760.0 / data.at("pressure").get<double>(),
            .X_feed = data.at("X_feed").get<double>(),
            .t_ads_start = data.at("t_ads_start").get<double>(),
            .t_ads_end = data.at("t_ads_end").get<double>(),
            .k_ads_smooth = data.at("k_ads_smooth").get<double>(),
            .dt =  dt_exp
        },
        .initial_guess = {
            .k_ads = data.at("initial_guess").at("k_ads").get<double>(),
            .k_des = data.at("initial_guess").at("k_des").get<double>(),
            .k_rxn = data.at("initial_guess").at("k_rxn").get<double>(),
            .S_tot = data.at("initial_guess").at("S_tot").get<double>(),
            .P_tot = data.at("initial_guess").at("Y_tot").get<double>()
        },
        .t_exp = t_exp,
        .X_exp = X_exp,
        .N_reactors = data.at("N_reactors").get<long>()
    };

    return input_data;
}

std::array<std::vector<double>, 10> solve_model(
    Model::FixedParameters const & fixed_parameters,
    Model::FittedParameters const & fitted_parameters,
    const long N_t, const long N_reactors) {

    Model model(fixed_parameters, fitted_parameters, N_reactors);

    std::array<std::vector<double>, 10> X;
    std::fill(X.begin(), X.end(), std::vector<double>(N_t));

    for (long i = 0; i < N_t; i ++) {

        model.do_step();

        for (int j = 0; j < 4; j ++)
            X[j][i] = model.get_Y()[model.xbase(N_reactors-1, j)];
        for (int j = 0; j < 5; j ++)
            X[j+4][i] = model.get_Y()[model.sbase(N_reactors-1, j)];

        X[9][i] = model.f_of_t(model.get_t());
    }

    return X;
}

int main(int argc, char ** argv) {

    if (argc < 2) {
        fmt::println(stderr, "Input file path must be provided as an argument");
        return EXIT_FAILURE;
    }

    std::filesystem::path input_file_path(argv[1]);
    std::filesystem::path output_file_path = input_file_path.parent_path() / "fitted.json";

    auto input_data = load_input_data(input_file_path);

    auto * cost = new ResidualFunctor(input_data.t0_exp, input_data.fixed_parameters, input_data.X_exp, input_data.N_reactors);

    // Initial guesses
    double theta[] = {
        S_LOG(input_data.initial_guess.k_ads),     // k_ads
        S_LOG(input_data.initial_guess.k_des),     // k_des
        S_LOG(input_data.initial_guess.k_rxn),     // k_rxn
        S_LOG(input_data.initial_guess.S_tot),     // S0
        S_LOG(input_data.initial_guess.P_tot)      // Y0
    };

    Model::FittedParameters fitted_parameters_0 {
        S_EXP(theta[0]), S_EXP(theta[1]), S_EXP(theta[2]), S_EXP(theta[3]), S_EXP(theta[4])
    };

    auto X0 = solve_model(input_data.fixed_parameters, fitted_parameters_0, input_data.t_exp.size(), input_data.N_reactors);

    ceres::Problem problem;
    problem.AddResidualBlock(cost, nullptr, theta);

    // Add constraints
#ifndef ENABLE_SCALING
    problem.SetParameterLowerBound(theta, 0, 0.0);
    problem.SetParameterLowerBound(theta, 1, 0.0);
    problem.SetParameterLowerBound(theta, 2, 0.0);
    problem.SetParameterLowerBound(theta, 3, 0.0);
    problem.SetParameterLowerBound(theta, 4, 0.0);
#endif //ENABLE_SCALING

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 200;  // increase from default 50
    // If your EvaluateAll is not thread-safe (CVODE usually isn't), force 1 thread:
    // options.gradient_check_numeric_derivative_relative_step_size = 1e-1;
    // options.check_gradients = true;
    options.num_threads = 1;

    // Tight tolerances - ONLY FOR FITTING ARTIFICIAL CURVES
    // options.function_tolerance   = 1e-12;   // cost change tolerance
    // options.gradient_tolerance   = 1e-14;   // gradient norm tolerance
    // options.parameter_tolerance  = 0.0;   // parameter step tolerance
    // options.use_nonmonotonic_steps = true;

    double cost_before;
    problem.Evaluate(ceres::Problem::EvaluateOptions(),
                     &cost_before, nullptr, nullptr, nullptr);
    std::cout << "Initial cost: " << cost_before << std::endl;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.FullReport() << "\n";

    Model::FittedParameters fitted_parameters {
        S_EXP(theta[0]), S_EXP(theta[1]), S_EXP(theta[2]), S_EXP(theta[3]), S_EXP(theta[4])
    };

    auto X = solve_model(input_data.fixed_parameters, fitted_parameters, input_data.t_exp.size(), input_data.N_reactors);

    // Calculate fit statistics
    double SSR = 0.0;
    for (long i = 0; i < X[0].size(); i ++) {
        double residual = (X[0][i] - input_data.X_exp[i]);
        SSR += residual * residual;
    }

    double sigma_sq = SSR / static_cast<double>(X[0].size() - 5);

    Eigen::Matrix2d J(X[0].size(), 5);

    for (long i = 0; i < X[0].size(); i ++) {
        for (long j = 0; j < 5; j ++) {
            J(i, j) = X[j+4][i];
        }
    }

    Eigen::Matrix2d cov = sigma_sq * (J.transpose() * J).inverse();
    for (long j = 0; j < 5; j ++) {
        fmt::println("Uncertainty: {}", sqrt(cov(j, j)));
    }

    // J_ij = df(x_i)/dk_j
    // cov = sigma_sq * (J^T*J)^-1
    // std. err j = sqrt(cov_jj)

    ordered_json out_data;
    out_data["solution"]["k_ads"] = S_EXP(theta[0]);
    out_data["solution"]["k_des"] = S_EXP(theta[1]);
    out_data["solution"]["k_rxn"] = S_EXP(theta[2]);
    out_data["solution"]["S_tot"] = S_EXP(theta[3]);
    out_data["solution"]["Y_tot"] = S_EXP(theta[4]);
    out_data["fitted_data"]["X0"] = X0[0];
    out_data["fitted_data"]["X"] = X[0];
    out_data["fitted_data"]["Xgs"] = X[1];
    out_data["fitted_data"]["Xs"] = X[2];
    out_data["fitted_data"]["P"] = X[3];

    // Outputting sensitivities...
    out_data["sensitivities"]["dXdk_ads"] = X[4];
    out_data["sensitivities"]["dXdk_des"] = X[5];
    out_data["sensitivities"]["dXdk_rxn"] = X[6];
    out_data["sensitivities"]["dXdS_tot"] = X[7];
    out_data["sensitivities"]["dXdY_tot"] = X[8];

    std::vector<double> uptake_rate(input_data.t_exp.size());

    for (int i = 0; i < uptake_rate.size(); i ++) {
        uptake_rate[i] = theta[0] * X[9][i] * X[1][i] * (theta[3] - X[2][i]) - input_data.fixed_parameters.R / 2.0 * theta[1] * X[2][i];
    }

    out_data["fitted_data"]["uptake_rate"] = uptake_rate;

    std::ofstream ofs(output_file_path);

    if (!ofs.good()) {
        fmt::println(stderr, "Error writing output to {}", output_file_path.string());
        exit(EXIT_FAILURE);
    }

    ofs << out_data.dump(4);

    return 0;
}