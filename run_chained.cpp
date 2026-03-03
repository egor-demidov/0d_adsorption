#include <vector>
#include <fstream>
#include <array>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <nlohmann/json.hpp>

#include "model_chained.h"

using json = nlohmann::json;
using ordered_json = nlohmann::ordered_json;

struct InputData {
    Model::FixedParameters fixed_parameters;
    Model::FittedParameters fitted_parameters;
    double t_tot;
    long N_t;
    long N_reactors;
};

InputData load_input_data(std::filesystem::path const & input_file_path) {

    std::ifstream ifs(input_file_path);

    if (!ifs.good()) {
        fmt::println(stderr, "Error reading input data from {}", input_file_path.string());
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

    long N_t = require_long(data, "N_t");
    long N_reactors = require_long(data, "N_reactors");

    auto fitted_parameters = require_object(data, "fitted_parameters");

    InputData input_data{
        .fixed_parameters = {
            .Di = require_double(data, "Di"),
            .R =  require_double(data, "R"),
            .L =  require_double(data, "L"),
            .F =  require_double(data, "F") * 760.0 /  require_double(data, "pressure"),
            .X_feed =  require_double(data, "X_feed"),
            .t_ads_start =  require_double(data, "t_ads_start"),
            .t_ads_end =  require_double(data, "t_ads_end"),
            .k_ads_smooth =  require_double(data, "k_ads_smooth"),
            .dt =  require_double(data, "t_tot") / static_cast<double>(N_t)
        },
        .fitted_parameters = {
            .k_ads = require_double(fitted_parameters, "k_ads"),
            .k_des = require_double(fitted_parameters, "k_des"),
            .k_rxn = require_double(fitted_parameters, "k_rxn"),
            .S_tot = require_double(fitted_parameters, "S_tot"),
            .P_tot = require_double(fitted_parameters, "Y_tot")
        },
        .N_t = N_t,
        .N_reactors = N_reactors
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
    std::filesystem::path output_file_path = input_file_path.parent_path() / "run.json";

    InputData input_data;

    try {
        input_data = load_input_data(input_file_path);
    } catch (std::runtime_error const & e) {
        fmt::println(stderr, "Error while parsing parameters: {}", e.what());
        return EXIT_FAILURE;
    }

    // Initial guesses
    Model::FittedParameters fitted_parameters {
            input_data.fitted_parameters.k_ads,     // k_ads
            input_data.fitted_parameters.k_des,     // k_des
            input_data.fitted_parameters.k_rxn,     // k_rxn
            input_data.fitted_parameters.S_tot,     // S0
            input_data.fitted_parameters.P_tot      // Y0
    };

    std::vector<double> ts(input_data.N_t);
    for (long i = 0; i < input_data.N_t; i ++) {
        ts[i] = static_cast<double>(i+1) * input_data.fixed_parameters.dt;
    }

    auto X = solve_model(input_data.fixed_parameters, fitted_parameters, input_data.N_t, input_data.N_reactors);

    // double h = 1.0e-16;
    // fitted_parameters.k_ads += h;

    // auto X_fwd = solve_model(input_data.fixed_parameters, fitted_parameters, input_data.N_t, input_data.N_reactors);

    // std::vector<double> dXfwd(X[0].size());

    // for (int i = 0; i < X[0].size(); i ++) {
        // dXfwd[i] = (X_fwd[0][i] - X[0][i]) / h;
    // }

    ordered_json out_data;
    out_data["solution"]["X"] = X[0];
    out_data["solution"]["t"] = ts;
    // out_data["fitted_data"]["X0"] = X0[0];
    // out_data["fitted_data"]["X"] = X[0];
    //
    // // Outputting sensitivities...
    out_data["sensitivities"]["dXdk_ads"] = X[4];
    out_data["sensitivities"]["dXdk_des"] = X[5];
    out_data["sensitivities"]["dXdk_rxn"] = X[6];
    out_data["sensitivities"]["dXdS_tot"] = X[7];
    out_data["sensitivities"]["dXdY_tot"] = X[8];
    // out_data["sensitivities"]["dXfwd"] = dXfwd;
    //
    // std::vector<double> uptake_rate(input_data.t_exp.size());
    //
    // for (int i = 0; i < uptake_rate.size(); i ++) {
    //     uptake_rate[i] = theta[0] * X[9][i] * X[1][i] * (theta[3] - X[2][i]) - input_data.fixed_parameters.R / 2.0 * theta[1] * X[2][i];
    // }
    //
    // out_data["fitted_data"]["uptake_rate"] = uptake_rate;
    //
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