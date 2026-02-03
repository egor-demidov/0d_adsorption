//
// Created by egor on 1/16/26.
//

#ifndef INC_0D_ADSORPTION_MODEL_H
#define INC_0D_ADSORPTION_MODEL_H

#include <sundials/sundials_types.h>     // sunrealtype
#include <sundials/sundials_context.h>   // SUNContext_Create, SUNContext_Free

#include <cvode/cvode.h>                 // CVode*
#include <nvector/nvector_serial.h>      // N_Vector, N_VNew_Serial
#include <sunmatrix/sunmatrix_dense.h>   // SUNMatrix_Dense
#include <sunlinsol/sunlinsol_dense.h>   // SUNLinSol_Dense
#include <sundials/sundials_math.h>      // SUNRabs

class Model {
public:
    struct FixedParameters {
        const sunrealtype Di;
        const sunrealtype R;
        const sunrealtype L;
        const sunrealtype F;
        const sunrealtype X_feed;
        const sunrealtype t_ads_start;
        const sunrealtype t_ads_end;
        const sunrealtype k_ads_smooth;
        const sunrealtype dt;
    };

    struct FittedParameters {
        const sunrealtype k_ads;
        const sunrealtype k_des;
        const sunrealtype k_rxn;
        const sunrealtype S_tot;
        const sunrealtype P_tot;
    };

    Model(FixedParameters const & fixed_parameters, FittedParameters const & fitted_parameters, int N_reactors);
    ~Model();

    [[nodiscard]]
    sunrealtype f_of_t(sunrealtype t) const {
        return 0.5 * (tanh((t - fixed_parameters_.t_ads_start) / fixed_parameters_.k_ads_smooth) - tanh((t - fixed_parameters_.t_ads_end) / fixed_parameters_.k_ads_smooth));
    }

    // n - index of reactor
    // i - index of parameter
    // returns index of sensitivity of X_n with respect to parameter n
    int sbase(int n, int j) const {
        return n * NX * (1 + NP) + NX * (1 + j);
    }

    // n - index of reactor
    // i - index of state
    int xbase(int n, int i) const {
        return n * NX * (1 + NP) + i;
    }

    void do_step();

    [[nodiscard]]
    const sunrealtype * get_Y() const {
        return N_VGetArrayPointer(y_);
    }

    [[nodiscard]]
    sunrealtype get_t() const {
        return t_out;
    }

private:
    struct DerivedParameters {
        const sunrealtype V;
        const sunrealtype k_diff;
        const sunrealtype a;
    };

    // ODE right hand side
    int rhs(sunrealtype t, N_Vector y, N_Vector ydot) const;
    // ODE Jacobial
    int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) const;

    static int rhs_cvode(sunrealtype t,
                       N_Vector y,
                       N_Vector ydot,
                       void* user_data) {
        // Recover the object
        auto* self = static_cast<Model*>(user_data);
        return self->rhs(t, y, ydot);
    }

    static int jac_cvode(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        // Recover the object
        auto* self = static_cast<Model*>(user_data);
        return self->jac(t, y, fy, J, tmp1, tmp2, tmp3);
    }

    // Parameter count
    const int N_reactors_;
    const int NP, NX, N;

    const FixedParameters fixed_parameters_;
    const FittedParameters fitted_parameters_;
    const DerivedParameters derived_parameters_;

    SUNContext sunctx_;
    N_Vector y_;
    N_Vector atol_;
    void* cvode_mem_;
    SUNMatrix A_;
    SUNLinearSolver LS_;

    sunrealtype t_out = 0.0;
};

#endif //INC_0D_ADSORPTION_MODEL_H