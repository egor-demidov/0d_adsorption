#include <vector>
#include <fstream>
#include <array>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <ceres/ceres.h>

#include <nlohmann/json.hpp>

// New JSON library used to control format of numbers
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>

#include "model_chained.h"

// #define ENABLE_SCALING

#ifdef ENABLE_SCALING
#define S_EXP(__X__) exp(__X__)
#define S_LOG(__X__) log(__X__)
#else //ENABLE_SCALING
#define S_EXP(__X__) (__X__)
#define S_LOG(__X__) (__X__)
#endif //ENABLE_SCALING

#ifdef EMSCRIPTEN
#include <emscripten/bind.h>
#endif //EMSCRIPTEN

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
            // fmt::println("{} {}", t0_exp_, model.get_t());
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

    auto require_object = [](const nlohmann::json& j, const std::string& key)
    -> const nlohmann::json& {

        if (!j.contains(key))
            throw std::runtime_error("Missing required object: '" + key + "'");

        if (!j.at(key).is_object())
            throw std::runtime_error("'" + key + "' must be an object");

        return j.at(key);
    };

    auto require_double = [](const nlohmann::json& j, const std::string& key) -> double {
        if (!j.contains(key)) {
            throw std::runtime_error("Missing required parameter: '" + key + "'");
        }
        if (!j.at(key).is_number()) {
            throw std::runtime_error("Parameter '" + key + "' must be a number");
        }
        return j.at(key).get<double>();
    };

    auto require_double_vector = [](const nlohmann::json& j, const std::string& key) -> std::vector<double> {
        if (!j.contains(key)) {
            throw std::runtime_error("Missing required parameter: '" + key + "'");
        }
        if (!j.at(key).is_array()) {
            throw std::runtime_error("Parameter '" + key + "' must be an array of numbers");
        }
        return j.at(key).get<std::vector<double>>();
    };

    auto require_long = [](const nlohmann::json& j, const std::string& key) -> long {
        if (!j.contains(key)) {
            throw std::runtime_error("Missing required parameter: '" + key + "'");
        }
        if (!j.at(key).is_number_integer()) {
            throw std::runtime_error("Parameter '" + key + "' must be an integer");
        }
        return j.at(key).get<long>();
    };

    auto experimental_data = require_object(data, "experimental_data");

    auto t_exp = require_double_vector(experimental_data, "t_exp");

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

    auto X_exp = require_double_vector(experimental_data, "X_exp");

    // std::rotate(X_exp.rbegin(), X_exp.rbegin() + 1, X_exp.rend());
    // std::rotate(X_exp.begin(), X_exp.begin() + 1, X_exp.end());

    auto initial_guess = require_object(data, "initial_guess");

    double pressure = require_double(data, "pressure");

    InputData input_data{
        .t0_exp = t0_exp,
        .fixed_parameters = {
            .Di = require_double(data, "Di") * 760.0 / pressure,
            .R =  require_double(data, "R"),
            .L =  require_double(data, "L"),
            .F =  require_double(data, "F") * 760.0 / pressure,
            .X_feed =  require_double(data, "X_feed"),
            .t_ads_start =  require_double(data, "t_ads_start"),
            .t_ads_end =  require_double(data, "t_ads_end"),
            .tau_sw_1 =  require_double(data, "tau_sw_1"),
            .tau_sw_2 =  require_double(data, "tau_sw_2"),
            .dt =  dt_exp
        },
        .initial_guess = {
            .k_ads = require_double(initial_guess, "k_ads"),
            .k_des = require_double(initial_guess, "k_des"),
            .k_rxn = require_double(initial_guess, "k_rxn"),
            .S_tot = require_double(initial_guess, "S_tot"),
            .P_tot = require_double(initial_guess, "Y_tot")
        },
        .t_exp = t_exp,
        .X_exp = X_exp,
        .N_reactors = require_long(data, "N_reactors")
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

std::vector<double> fit_model(Model::FixedParameters fixed_parameters,
    const double t0_exp,
    std::vector<double> X_exp,
    const long N_reactors,
    std::vector<double> theta_vec) {

    auto * cost = new ResidualFunctor(t0_exp, fixed_parameters, X_exp, N_reactors);

    double * theta = &theta_vec[0];

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

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.FullReport() << "\n";

    return theta_vec;
}

using PairDD = std::pair<double, double>;

struct FitResult {
    std::vector<double> X_sol, X_sol_0, X_exp, t_exp;
    std::pair<double, double> k_ads;
    std::pair<double, double> k_des;
    std::pair<double, double> k_rxn;
    std::pair<double, double> S_tot;
    std::pair<double, double> Y_tot;
};

FitResult run_fitting_web_app(Model::FixedParameters fixed_parameters,
    const double t0_exp,
    long t_exp_size,
    std::vector<double> X_exp,
    const long N_reactors,
    std::vector<double> theta_vec) {

    Model::FittedParameters fitted_parameters_0 {
        S_EXP(theta_vec[0]), S_EXP(theta_vec[1]), S_EXP(theta_vec[2]), S_EXP(theta_vec[3]), S_EXP(theta_vec[4])
    };

    auto X0 = solve_model(fixed_parameters, fitted_parameters_0, t_exp_size, N_reactors);

    // Perform the fitting
    theta_vec = fit_model(fixed_parameters, t0_exp, X_exp, N_reactors, theta_vec);

    Model::FittedParameters fitted_parameters {
        S_EXP(theta_vec[0]), S_EXP(theta_vec[1]), S_EXP(theta_vec[2]), S_EXP(theta_vec[3]), S_EXP(theta_vec[4])
    };

    auto X = solve_model(fixed_parameters, fitted_parameters, t_exp_size, N_reactors);

    // Calculate fit statistics
    double SSR = 0.0;
    for (long i = 0; i < X[0].size(); i ++) {
        double residual = (X[0][i] - X_exp[i]);
        SSR += residual * residual;
    }

    double sigma_sq = SSR / static_cast<double>(X[0].size() - 5);

    Eigen::MatrixXd J(X[0].size(), 5);

    for (long i = 0; i < X[0].size(); i ++) {
        for (long j = 0; j < 5; j ++) {
            J(i, j) = X[j+4][i];
        }
    }

    Eigen::MatrixXd cov = sigma_sq * (J.transpose() * J).inverse();

    std::vector<double> t_exp(X_exp.size());
    for (long i = 0; i < t_exp.size(); i ++) {
        t_exp[i] = t0_exp + static_cast<double>(i) * fixed_parameters.dt;
    }

    return {
        .X_sol = X[0],
        .X_sol_0 = X0[0],
        .X_exp = X_exp,
        .t_exp = t_exp,
        .k_ads = std::make_pair(S_EXP(theta_vec[0]), S_EXP(sqrt(cov(0, 0)))),
        .k_des = std::make_pair(S_EXP(theta_vec[1]), S_EXP(sqrt(cov(1, 1)))),
        .k_rxn = std::make_pair(S_EXP(theta_vec[2]), S_EXP(sqrt(cov(2, 2)))),
        .S_tot = std::make_pair(S_EXP(theta_vec[3]), S_EXP(sqrt(cov(3, 3)))),
        .Y_tot = std::make_pair(S_EXP(theta_vec[4]), S_EXP(sqrt(cov(4, 4))))
    };
}

#ifdef EMSCRIPTEN
EMSCRIPTEN_BINDINGS(0d_adsorption_fit_chained) {
    emscripten::register_vector<double>("VectorDouble");
    emscripten::function("run_fitting_web_app", &run_fitting_web_app);

    emscripten::value_object<Model::FixedParameters>("FixedParameters")
        .field("Di", &Model::FixedParameters::Di)
        .field("R", &Model::FixedParameters::R)
        .field("L", &Model::FixedParameters::L)
        .field("F", &Model::FixedParameters::F)
        .field("X_feed", &Model::FixedParameters::X_feed)
        .field("t_ads_start", &Model::FixedParameters::t_ads_start)
        .field("t_ads_end", &Model::FixedParameters::t_ads_end)
        .field("k_ads_smooth", &Model::FixedParameters::k_ads_smooth)
        .field("dt", &Model::FixedParameters::dt);

    emscripten::value_object<PairDD>("PairDD")
        .field("first", &PairDD::first)
        .field("second", &PairDD::second);

    emscripten::value_object<FitResult>("FitResult")
        .field("X_sol", &FitResult::X_sol)
        .field("X_sol_0", &FitResult::X_sol_0)
        .field("X_exp", &FitResult::X_exp)
        .field("t_exp", &FitResult::t_exp)
        .field("k_ads", &FitResult::k_ads)
        .field("k_des", &FitResult::k_des)
        .field("k_rxn", &FitResult::k_rxn)
        .field("S_tot", &FitResult::S_tot)
        .field("Y_tot", &FitResult::Y_tot);

    // emscripten::value_object<Model::FittedParameters>("FittedParameters")
    //     .field("k_ads", &Model::FittedParameters::k_ads)
    //     .field("k_des", &Model::FittedParameters::k_des)
    //     .field("k_rxn", &Model::FittedParameters::k_rxn)
    //     .field("S_tot", &Model::FittedParameters::S_tot)
    //     .field("P_tot", &Model::FittedParameters::P_tot);
    //
    // emscripten::value_object<InputData>("InputData")
    //     .field("t0_exp", &InputData::t0_exp)
    //     .field("fixed_parameters", &InputData::fixed_parameters)
    //     .field("initial_guess", &InputData::initial_guess)
    //     .field("t_exp", &InputData::t_exp)
    //     .field("X_exp", &InputData::X_exp)
    //     .field("N_reactors", &InputData::N_reactors);


    // emscripten::class_<InputData>("InputData")
    //     .constructor<std::string, double, double, double>()
    //     .property("name_readonly", &Component::get_name)
    //     .property("density_readonly", &Component::get_density)
    //     .property("molecular_weight_readonly", &Component::get_molecular_weight)
    //     .function("p_sat", &Component::get_p_sat)
    //     .property("molar_volume_readonly", &Component::get_molar_volume);
    //
    // emscripten::value_object<SingleComponentCapillaryCondensationRun::Solution>("CondensationSolution")
    //     .field("time", &SingleComponentCapillaryCondensationRun::Solution::time)
    //     .field("condensate_volume", &SingleComponentCapillaryCondensationRun::Solution::condensate_volume)
    //     .field("condensate_volume_fraction", &SingleComponentCapillaryCondensationRun::Solution::condensate_volume_fraction)
    //     .field("uniform_to_capillary_ratio", &SingleComponentCapillaryCondensationRun::Solution::uniform_to_capillary_ratio)
    //     .field("capillary_filling_angle", &SingleComponentCapillaryCondensationRun::Solution::capillary_filling_angle)
    //     .field("uniform_coating_thickness", &SingleComponentCapillaryCondensationRun::Solution::uniform_coating_thickness)
    //     .field("capillary_condensate_volume", &SingleComponentCapillaryCondensationRun::Solution::capillary_condensate_volume)
    //     .field("uniform_condensate_volume", &SingleComponentCapillaryCondensationRun::Solution::uniform_condensate_volume);
    //
    // emscripten::value_object<CapillaryCondensationResult>("CapillaryCondensationResult")
    //     .field("solution", &CapillaryCondensationResult::solution)
    //     .field("n_steps", &CapillaryCondensationResult::n_steps)
    //     .field("dt", &CapillaryCondensationResult::dt)
    //     .field("chi", &CapillaryCondensationResult::chi)
    //     .field("capillary_condensation_threshold", &CapillaryCondensationResult::capillary_condensation_threshold);
}
#endif //EMSCRIPTEN

#ifndef EMSCRIPTEN
int main(int argc, char ** argv) {

    if (argc < 2) {
        fmt::println(stderr, "Input file path must be provided as an argument");
        return EXIT_FAILURE;
    }

    std::filesystem::path input_file_path(argv[1]);
    std::filesystem::path output_file_path = input_file_path.parent_path() / "fitted.json";

    InputData input_data;

    try {
        input_data = load_input_data(input_file_path);
    } catch (std::runtime_error const & e) {
        fmt::println(stderr, "Error while parsing parameters: {}", e.what());
        return EXIT_FAILURE;
    }

    // Initial guesses
    std::vector<double> theta = {
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

    // Perform the fitting
    theta = fit_model(input_data.fixed_parameters, input_data.t0_exp, input_data.X_exp, input_data.N_reactors, theta);

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

    Eigen::MatrixXd J(X[0].size(), 5);

    for (long i = 0; i < X[0].size(); i ++) {
        for (long j = 0; j < 5; j ++) {
            J(i, j) = X[j+4][i];
        }
    }

    Eigen::MatrixXd cov = sigma_sq * (J.transpose() * J).inverse();

    // J_ij = df(x_i)/dk_j
    // cov = sigma_sq * (J^T*J)^-1
    // std. err j = sqrt(cov_jj)

    std::vector<double> uptake_rate(input_data.t_exp.size());

    for (int i = 0; i < uptake_rate.size(); i ++) {
        uptake_rate[i] = theta[0] * X[9][i] * X[1][i] * (theta[3] - X[2][i]) - theta[1] * X[2][i];
    }

    rapidjson::StringBuffer sb;
    rapidjson::PrettyWriter<rapidjson::StringBuffer> w(sb);
    w.SetIndent(' ', 4);

    w.StartObject(); // root
    w.Key("solution");
    w.StartObject(); // solution

    w.Key("k_ads");
    auto k_ads = fmt::format("{:.06e}", S_EXP(theta[0]));
    w.RawValue(k_ads.c_str(), static_cast<rapidjson::SizeType>(k_ads.size()), rapidjson::kNumberType);

    w.Key("k_des");
    auto k_des = fmt::format("{:.06e}", S_EXP(theta[1]));
    w.RawValue(k_des.c_str(), static_cast<rapidjson::SizeType>(k_des.size()), rapidjson::kNumberType);

    w.Key("k_rxn");
    auto k_rxn = fmt::format("{:.06e}", S_EXP(theta[2]));
    w.RawValue(k_rxn.c_str(), static_cast<rapidjson::SizeType>(k_rxn.size()), rapidjson::kNumberType);

    w.Key("S_tot");
    auto S_tot = fmt::format("{:.06e}", S_EXP(theta[3]));
    w.RawValue(S_tot.c_str(), static_cast<rapidjson::SizeType>(S_tot.size()), rapidjson::kNumberType);

    w.Key("Y_tot");
    auto Y_tot = fmt::format("{:.06e}", S_EXP(theta[4]));
    w.RawValue(Y_tot.c_str(), static_cast<rapidjson::SizeType>(Y_tot.size()), rapidjson::kNumberType);

    w.EndObject(); // solution

    w.Key("standard_error");
    w.StartObject(); // standard_error

    w.Key("k_ads");
    k_ads = fmt::format("{:.06e}", S_EXP(sqrt(cov(0, 0))));
    w.RawValue(k_ads.c_str(), static_cast<rapidjson::SizeType>(k_ads.size()), rapidjson::kNumberType);

    w.Key("k_des");
    k_des = fmt::format("{:.06e}", S_EXP(sqrt(cov(1, 1))));
    w.RawValue(k_des.c_str(), static_cast<rapidjson::SizeType>(k_des.size()), rapidjson::kNumberType);

    w.Key("k_rxn");
    k_rxn = fmt::format("{:.06e}", S_EXP(sqrt(cov(2, 2))));
    w.RawValue(k_rxn.c_str(), static_cast<rapidjson::SizeType>(k_rxn.size()), rapidjson::kNumberType);

    w.Key("S_tot");
    S_tot = fmt::format("{:.06e}", S_EXP(sqrt(cov(3, 3))));
    w.RawValue(S_tot.c_str(), static_cast<rapidjson::SizeType>(S_tot.size()), rapidjson::kNumberType);

    w.Key("Y_tot");
    Y_tot = fmt::format("{:.06e}", S_EXP(sqrt(cov(4, 4))));
    w.RawValue(Y_tot.c_str(), static_cast<rapidjson::SizeType>(Y_tot.size()), rapidjson::kNumberType);

    w.EndObject(); // standard_error

    w.Key("fitted_data");
    w.StartObject(); // fitted_data

    w.Key("X0");
    w.StartArray(); // X0
    for (auto v : X0[0])
        w.Double(v);
    w.EndArray(); // X0

    w.Key("X");
    w.StartArray(); // X
    for (auto v : X[0])
        w.Double(v);
    w.EndArray(); // X

    w.Key("Xgs");
    w.StartArray(); // Xgs
    for (auto v : X[1])
        w.Double(v);
    w.EndArray(); // Xgs

    w.Key("Xs");
    w.StartArray(); // Xs
    for (auto v : X[2])
        w.Double(v);
    w.EndArray(); // Xs

    w.Key("P");
    w.StartArray(); // P
    for (auto v : X[3])
        w.Double(v);
    w.EndArray(); // P

    w.Key("uptake_rate");
    w.StartArray(); // uptake_rate
    for (auto v : uptake_rate)
        w.Double(v);
    w.EndArray(); // uptake_rate

    w.EndObject(); // fitted_data

    w.Key("sensitivities");
    w.StartObject(); // sensitivities

    w.Key("dXdk_ads");
    w.StartArray(); // dXdk_ads
    for (auto v :  X[4])
        w.Double(v);
    w.EndArray(); // dXdk_ads

    w.Key("dXdk_des");
    w.StartArray(); // dXdk_des
    for (auto v :  X[5])
        w.Double(v);
    w.EndArray(); // dXdk_des

    w.Key("dXdk_rxn");
    w.StartArray(); // dXdk_rxn
    for (auto v :  X[6])
        w.Double(v);
    w.EndArray(); // dXdk_rxn

    w.Key("dXdS_tot");
    w.StartArray(); // dXdS_tot
    for (auto v :  X[7])
        w.Double(v);
    w.EndArray(); // dXdS_tot

    w.Key("dXdY_tot");
    w.StartArray(); // dXdY_tot
    for (auto v :  X[8])
        w.Double(v);
    w.EndArray(); // dXdY_tot

    w.EndObject(); // sensitivities

    w.EndObject(); // root

    std::ofstream ofs(output_file_path);

    if (!ofs.good()) {
        fmt::println(stderr, "Error writing output to {}", output_file_path.string());
        exit(EXIT_FAILURE);
    }

    ofs << sb.GetString();

    return 0;
}
#endif //EMSCRIPTEN
