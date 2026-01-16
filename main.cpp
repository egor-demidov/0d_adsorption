#include <vector>
#include <fstream>
#include <array>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include "model.h"

int main() {
    Model::FixedParameters fixed_params{
        .Di = 40.0,
        .R = 1.56 / 2.0,
        .L = 2.0,
        .F = 2.493333333333333 * 760.0 / 1.965,
        .X_in = 29424160260.695,
        .t_ads_start =  266.3915963445649,
        .t_ads_end = 543.5885793335225,
        .k_ads_smooth = 2.645834006658198,
        .dt = 1.0
    };

    Model::FittedParameters fitted_params{
        .k_ads  = 3.827187367718609e-12,
        .k_des  = 0.04253902879091293,
        .k_rxn  = 2.0493834170301142e-16,
        .S_tot  = 34482079719017.1,
        .P_tot  = 88448995701717.02
    };

    Model model(fixed_params, fitted_params);

    const long N_t = 800;
    std::vector<double> ts(N_t), Xs(N_t);
    std::vector<std::array<double, 5>> dXs(N_t, std::array<double, 5>{});

    ts[0] = 0.0;
    Xs[0] = fixed_params.X_in;
    for (long i = 1; i < N_t; i ++) {

        model.do_step();

        ts[i] = model.get_t();
        Xs[i] = model.get_Y()[0];
        for (int j = 0; j < 5; j ++)
            dXs[i][j] = model.get_Y()[Model::sbase(j)];
    }

    std::ofstream ofs("out_data.csv");

    fmt::println(ofs, "t, X, dX/dk_ads, dX/dk_des, dX/dk_rxn, dX/dS0, dX/dY0");

    double multipliers[] = {fitted_params.k_ads, fitted_params.k_des, fitted_params.k_rxn, fitted_params.S_tot, fitted_params.P_tot};

    for (long i = 0; i < N_t; i ++) {
        fmt::print(ofs, "{}, {}", ts[i], Xs[i]);
        for (int j = 0; j < 5; j ++) {
            fmt::print(ofs, ", {}", dXs[i][j]*multipliers[j]);
        }
        fmt::println(ofs, "");
    }

    return 0;
}