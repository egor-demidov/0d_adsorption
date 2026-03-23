// cvode_augmented.cpp
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <fstream>

#include <sundials/sundials_types.h>     // sunrealtype
#include <sundials/sundials_context.h>   // SUNContext_Create, SUNContext_Free

#include <cvode/cvode.h>                 // CVode*
#include <nvector/nvector_serial.h>      // N_Vector, N_VNew_Serial
#include <sunmatrix/sunmatrix_dense.h>   // SUNMatrix_Dense
#include <sunlinsol/sunlinsol_dense.h>   // SUNLinSol_Dense
#include <sundials/sundials_math.h>      // SUNRabs

// Banded matrix definitions
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_band.h>

#include <fmt/format.h>

#include "model_chained.h"

Model::Model(FixedParameters const & fixed_parameters, FittedParameters const & fitted_parameters, int N_reactors)
  : N_reactors_(N_reactors)
  , NP(5)
  , NX(4)
  , N(N_reactors * NX * (1 + NP))
  , B(NX * (1 + NP))
  , mu(B - 1)
  , ml(2*B - 1)
  , fixed_parameters_(fixed_parameters)
  , fitted_parameters_(fitted_parameters)
  , derived_parameters_{
    .V = M_PI * fixed_parameters.R * fixed_parameters.R * fixed_parameters.L,
    .k_diff = 3.66 * fixed_parameters.Di / fixed_parameters.R / fixed_parameters.R,
    .a = fixed_parameters.F / (M_PI * fixed_parameters.R * fixed_parameters.R * fixed_parameters.L / static_cast<double>(N_reactors))
  }
  , sunctx_(nullptr)
  , y_(nullptr)
  , atol_(nullptr)
  , cvode_mem_(nullptr)
  , A_(nullptr)
  , LS_(nullptr)
{
  // Initial conditions for x:
  sunrealtype X0=fixed_parameters.X_feed, Xgs0=fixed_parameters.X_feed, Xs0=0, P0=0;

  // Create a sundials context
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx_) != 0) {
    std::fprintf(stderr, "SUNContext_Create failed\n");
    exit(EXIT_FAILURE);
  }

  // Create N_Vector y0 and set the initial conditions
  y_ = N_VNew_Serial(N, sunctx_);
  sunrealtype* Y = N_VGetArrayPointer(y_);
  for (int i=0;i<N;++i) Y[i]=0.0;

  // Set the initial conditions for each reactor
  for (int n = 0; n < N_reactors_; n ++) {
    Y[xbase(n, 0)] = X0;
    Y[xbase(n, 1)] = Xgs0;
    Y[xbase(n, 2)] = Xs0;
    Y[xbase(n, 3)] = P0;
  }

  // Sensitivities s(0)=0 already set.

  // Create CVODE memory
  cvode_mem_ = CVodeCreate(CV_BDF, sunctx_);  // stiff -> BDF
  if (!cvode_mem_) { std::fprintf(stderr,"CVodeCreate failed\n"); exit(EXIT_FAILURE); }

  int flag = CVodeInit(cvode_mem_, rhs_cvode, t_out, y_);
  if (flag != CV_SUCCESS) { std::fprintf(stderr,"CVodeInit failed\n"); exit(EXIT_FAILURE); }

  CVodeSetMaxNumSteps(cvode_mem_, 200000); // or 1000000 if needed
  CVodeSetMaxStep(cvode_mem_, fixed_parameters_.dt);
  CVodeSetInitStep(cvode_mem_, fixed_parameters_.dt);

  flag = CVodeSetUserData(cvode_mem_, this);
  if (flag != CV_SUCCESS) { std::fprintf(stderr,"CVodeSetUserData failed\n"); exit(EXIT_FAILURE); }

  // tolerances
  sunrealtype reltol = 1e-8;

  atol_ = N_VNew_Serial(N, sunctx_);
  auto* a = N_VGetArrayPointer(atol_);

  // Set absolute tolerances
  // Sensitivities: MUCH looser at first
  for (int i = 0; i < N; ++i) a[i] = 1e6;

  // States: scale to your magnitude (~1e10 for X, Xgs)
  for (int n = 0; n < N_reactors_; n ++) {
    a[xbase(n, 0)] = 1e2;   // X
    a[xbase(n, 1)] = 1e2;   // Xgs
    a[xbase(n, 2)] = 1e-6;  // Xs (if small)
    a[xbase(n, 3)] = 1e-6;  // P  (if small)
  }

  flag = CVodeSVtolerances(cvode_mem_, reltol, atol_);
  if (flag != CV_SUCCESS) { std::fprintf(stderr,"CVodeSStolerances failed\n"); exit(EXIT_FAILURE); }

  // Dense linear solver + dense Jacobian
  // A_ = SUNDenseMatrix(N, N, sunctx_);;
  // LS_ = SUNLinSol_Dense(y_, A_, sunctx_);

  A_  = SUNBandMatrix(N, mu, ml, sunctx_);
  LS_ = SUNLinSol_Band(y_, A_, sunctx_);

  flag = CVodeSetLinearSolver(cvode_mem_, LS_, A_);
  if (flag != CV_SUCCESS) { std::fprintf(stderr,"CVodeSetLinearSolver failed\n"); exit(EXIT_FAILURE); }


  flag = CVodeSetJacFn(cvode_mem_, jac_cvode);
  // flag = CVodeSetJacFn(cvode_mem_, nullptr);
  if (flag != CV_SUCCESS) { std::fprintf(stderr,"CVodeSetJacFn failed\n"); exit(EXIT_FAILURE); }
}

Model::~Model() {
  // Cleanup
  CVodeFree(&cvode_mem_);
  SUNLinSolFree(LS_);
  SUNMatDestroy(A_);
  N_VDestroy(y_);
  N_VDestroy(atol_);

  SUNContext_Free(&sunctx_);
}

void Model::reset_model(FittedParameters const & new_fitted_parameters) {
  // Reset time
  t_out = 0.0;
  CVodeReInit(cvode_mem_, t_out, y_);

  // Overwrite the fitted parameters
  fitted_parameters_ = new_fitted_parameters;
}

int Model::rhs(sunrealtype t, N_Vector y, N_Vector ydot) const {
  const sunrealtype k_de = fitted_parameters_.k_des;
  const sunrealtype k_r  = fitted_parameters_.k_rxn;
  const sunrealtype S    = fitted_parameters_.S_tot;
  const sunrealtype PT   = fitted_parameters_.P_tot;

  const sunrealtype k_d  = derived_parameters_.k_diff;
  const sunrealtype a    = derived_parameters_.a;
  const sunrealtype Xfeed  = fixed_parameters_.X_feed;
  const sunrealtype R    = fixed_parameters_.R;
  const sunrealtype ft   = f_of_t(t);
  const sunrealtype k_a_eff = fitted_parameters_.k_ads * ft;
  // const sunrealtype k_a_eff = fitted_parameters_.k_ads;

  sunrealtype* Y  = N_VGetArrayPointer(y);
  sunrealtype* dY = N_VGetArrayPointer(ydot);

  for (int n = 0; n < N_reactors_; n ++) {

    const sunrealtype Xprev = (n == 0 ? Xfeed : Y[xbase(n-1, 0)]);

    // States
    const sunrealtype X   = Y[xbase(n, 0)];
    const sunrealtype Xgs = Y[xbase(n, 1)];
    const sunrealtype Xs  = Y[xbase(n, 2)];
    const sunrealtype P   = Y[xbase(n, 3)];

    // ---- f(x,theta) ----
    const sunrealtype f1 = a*(Xprev - X) - k_d*(X - Xgs);
    const sunrealtype f2 = k_d*(X - Xgs) - (2.0/R)*k_a_eff*Xgs*(S - Xs) + (2.0/R)*k_de*Xs;
    const sunrealtype f3 = k_a_eff*Xgs*(S - Xs) - k_de*Xs - k_r*Xs*(PT - P);
    const sunrealtype f4 = k_r*Xs*(PT - P);

    dY[xbase(n, 0)] = f1;
    dY[xbase(n, 1)] = f2;
    dY[xbase(n, 2)] = f3;
    dY[xbase(n, 3)] = f4;

    // ---- A = df/dx (4x4) ----
    const sunrealtype A11 = -a - k_d;
    const sunrealtype A12 =  k_d;

    const sunrealtype A21 =  k_d;
    const sunrealtype A22 = -k_d - (2.0/R)*k_a_eff*(S - Xs);
    const sunrealtype A23 =  (2.0/R)*k_a_eff*Xgs + (2.0/R)*k_de;

    const sunrealtype A32 =  k_a_eff*(S - Xs);
    const sunrealtype A33 = -k_a_eff*Xgs - k_de - k_r*(PT - P);
    const sunrealtype A34 =  k_r*Xs;

    const sunrealtype A43 =  k_r*(PT - P);
    const sunrealtype A44 = -k_r*Xs;

    // For each parameter, build b^(j) = df/dtheta_j and compute sdot = A*s + b
    for (int j = 0; j < NP; ++j) {
      const int b = sbase(n, j);
      const sunrealtype sX   = Y[b+0];
      const sunrealtype sXgs = Y[b+1];
      const sunrealtype sXs  = Y[b+2];
      const sunrealtype sP   = Y[b+3];

      // A*s
      sunrealtype As0 = A11*sX + A12*sXgs;

      // add upstream sensitivity coupling through Xprev = X_{n-1}
      if (n > 0) {
        const sunrealtype sX_prev = Y[sbase(n-1, j) + 0];
        As0 += a * sX_prev;
      }

      sunrealtype As1 = A21*sX + A22*sXgs + A23*sXs;
      sunrealtype As2 = A32*sXgs + A33*sXs + A34*sP;
      sunrealtype As3 = A43*sXs + A44*sP;

      // b^(j)
      sunrealtype b0=0, b1=0, b2=0, b3=0;

      switch (j) {
        case 0: // k_ads
          b1 = -(2.0/R)*Xgs*(S - Xs)*ft;
          b2 =  (1.0)*Xgs*(S - Xs)*ft;
          break;
        case 1: // k_des
          b1 =  (2.0/R)*Xs;
          b2 = -(1.0)*Xs;
          break;
        case 2: // k_rxn
          b2 = -Xs*(PT - P);
          b3 =  Xs*(PT - P);
          break;
        case 3: // S_tot
          b1 = -(2.0/R)*k_a_eff*Xgs;
          b2 =  k_a_eff*Xgs;
          break;
        case 4: // P_tot
          b2 = -k_r*Xs;
          b3 =  k_r*Xs;
          break;
        default:
          break;
      }

      dY[b+0] = As0 + b0;
      dY[b+1] = As1 + b1;
      dY[b+2] = As2 + b2;
      dY[b+3] = As3 + b3;
    }

  }

  return 0;

}

int Model::jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) const
{
  (void)fy; (void)tmp1; (void)tmp2; (void)tmp3;

  const sunrealtype k_de = fitted_parameters_.k_des;
  const sunrealtype k_r  = fitted_parameters_.k_rxn;
  const sunrealtype S    = fitted_parameters_.S_tot;
  const sunrealtype PT   = fitted_parameters_.P_tot;

  const sunrealtype k_d  = derived_parameters_.k_diff;
  const sunrealtype a    = derived_parameters_.a;
  const sunrealtype R    = fixed_parameters_.R;

  const sunrealtype ft     = f_of_t(t);
  const sunrealtype k_a_eff = fitted_parameters_.k_ads * ft;

  sunrealtype* Y = N_VGetArrayPointer(y);

  // auto setJ = [&](int row, int col, sunrealtype val) {
  //   SM_ELEMENT_D(J, row, col) = val;
  // };
  auto setJ = [&](sunindextype row, sunindextype col, sunrealtype val) {
    // Only set if (row,col) lies inside the stored band.
    const sunindextype r = row;
    const sunindextype c = col;
    if ((c <= r + mu) && (r <= c + ml)) {
      SM_ELEMENT_B(J, r, c) = val;
    }
  };

  // Zero J
  SUNMatZero(J);

  // Loop reactors
  for (int n = 0; n < N_reactors_; ++n) {

    // State indices for this reactor
    const int xX   = xbase(n, 0);
    const int xXgs = xbase(n, 1);
    const int xXs  = xbase(n, 2);
    const int xP   = xbase(n, 3);

    // Current reactor states
    const sunrealtype X   = Y[xX];
    const sunrealtype Xgs = Y[xXgs];
    const sunrealtype Xs  = Y[xXs];
    const sunrealtype P   = Y[xP];

    // ---- A(n) = df_n/dx_n ----
    const sunrealtype A11 = -a - k_d;
    const sunrealtype A12 =  k_d;

    const sunrealtype A21 =  k_d;
    const sunrealtype A22 = -k_d - (2.0/R)*k_a_eff*(S - Xs);
    const sunrealtype A23 =  (2.0/R)*k_a_eff*Xgs + (2.0/R)*k_de;

    const sunrealtype A32 =  k_a_eff*(S - Xs);
    const sunrealtype A33 = -k_a_eff*Xgs - k_de - k_r*(PT - P);
    const sunrealtype A34 =  k_r*Xs;

    const sunrealtype A43 =  k_r*(PT - P);
    const sunrealtype A44 = -k_r*Xs;

    // ---- State block (x_n wrt x_n) ----
    setJ(xX,   xX,   A11);  setJ(xX,   xXgs, A12);

    setJ(xXgs, xX,   A21);  setJ(xXgs, xXgs, A22);  setJ(xXgs, xXs, A23);

    setJ(xXs,  xXgs, A32);  setJ(xXs,  xXs,  A33);  setJ(xXs,  xP,  A34);

    setJ(xP,   xXs,  A43);  setJ(xP,   xP,   A44);

    // ---- Chain coupling in state equation: d(dX_n)/d(X_{n-1}) = a ----
    if (n > 0) {
      const int xX_prev = xbase(n-1, 0);
      setJ(xX, xX_prev, a);
    }

    // ---- Sensitivity blocks ----
    for (int j = 0; j < NP; ++j) {
      const int b = sbase(n, j);

      const sunrealtype sX   = Y[b+0];
      const sunrealtype sXgs = Y[b+1];
      const sunrealtype sXs  = Y[b+2];
      const sunrealtype sP   = Y[b+3];

      // Diagonal sensitivity block: d(sdot)/d(s) = A(n)
      setJ(b+0, b+0, A11); setJ(b+0, b+1, A12);

      setJ(b+1, b+0, A21); setJ(b+1, b+1, A22); setJ(b+1, b+2, A23);

      setJ(b+2, b+1, A32); setJ(b+2, b+2, A33); setJ(b+2, b+3, A34);

      setJ(b+3, b+2, A43); setJ(b+3, b+3, A44);

      // Chain coupling in sensitivity equation for sX:
      // sdotX_n = ... + a*sX_{n-1}
      if (n > 0) {
        setJ(b+0, sbase(n-1, j) + 0, a);
      }

      // ---- dA/d(state)*s terms (same as single reactor, local) ----
      // (dA/dXgs)*s = [0, (2/R)k_a*sXs, -k_a*sXs, 0]^T
      const sunrealtype dA_Xgs_0 = 0.0;
      const sunrealtype dA_Xgs_1 = (2.0/R)*k_a_eff*sXs;
      const sunrealtype dA_Xgs_2 = -k_a_eff*sXs;
      const sunrealtype dA_Xgs_3 = 0.0;

      // (dA/dXs)*s = [0, (2/R)k_a*sXgs, -k_a*sXgs, 0]^T
      const sunrealtype dA_Xs_0  = 0.0;
      const sunrealtype dA_Xs_1  = (2.0/R)*k_a_eff*sXgs;
      const sunrealtype dA_Xs_2  = -k_a_eff*sXgs;
      const sunrealtype dA_Xs_3  = 0.0;

      // (dA/dP)*s = [0, 0, k_r*sXs, -k_r*sXs]^T
      const sunrealtype dA_P_0 = 0.0;
      const sunrealtype dA_P_1 = 0.0;
      const sunrealtype dA_P_2 =  k_r*sXs;
      const sunrealtype dA_P_3 = -k_r*sXs;

      // ---- db/d(state) depends on parameter j ----
      sunrealtype db_Xgs_0=0, db_Xgs_1=0, db_Xgs_2=0, db_Xgs_3=0;
      sunrealtype db_Xs_0 =0, db_Xs_1 =0, db_Xs_2 =0, db_Xs_3 =0;
      sunrealtype db_P_0  =0, db_P_1  =0, db_P_2  =0, db_P_3  =0;

      switch (j) {
        case 0: // k_ads0
          db_Xgs_1 = -(2.0/R)*(S - Xs)*ft;
          db_Xgs_2 =  (S - Xs)*ft;
          db_Xs_1  =  (2.0/R)*Xgs*ft;
          db_Xs_2  = -Xgs*ft;
          break;
        case 1: // k_des
          db_Xs_1  =  (2.0/R);
          db_Xs_2  = -(1.0);
          break;
        case 2: // k_rxn
          db_Xs_2  = -(PT - P);
          db_Xs_3  =  (PT - P);
          db_P_2   =  Xs;
          db_P_3   = -Xs;
          break;
        case 3: // S_tot
          db_Xgs_1 = -(2.0/R)*k_a_eff;
          db_Xgs_2 =  k_a_eff;
          break;
        case 4: // P_tot
          db_Xs_2  = -k_r;
          db_Xs_3  =  k_r;
          break;
        default:
          break;
      }

      // ---- Fill D(n,j) into columns of THIS reactor's state block ----
      // Column xX is zero (A and b do not depend on X explicitly), so skip.

      // Column Xgs
      setJ(b+0, xXgs, dA_Xgs_0 + db_Xgs_0);
      setJ(b+1, xXgs, dA_Xgs_1 + db_Xgs_1);
      setJ(b+2, xXgs, dA_Xgs_2 + db_Xgs_2);
      setJ(b+3, xXgs, dA_Xgs_3 + db_Xgs_3);

      // Column Xs
      setJ(b+0, xXs,  dA_Xs_0  + db_Xs_0);
      setJ(b+1, xXs,  dA_Xs_1  + db_Xs_1);
      setJ(b+2, xXs,  dA_Xs_2  + db_Xs_2);
      setJ(b+3, xXs,  dA_Xs_3  + db_Xs_3);

      // Column P
      setJ(b+0, xP,   dA_P_0   + db_P_0);
      setJ(b+1, xP,   dA_P_1   + db_P_1);
      setJ(b+2, xP,   dA_P_2   + db_P_2);
      setJ(b+3, xP,   dA_P_3   + db_P_3);
    }
  }

  return 0;
}

void Model::do_step() {
  t_out += fixed_parameters_.dt;
  sunrealtype t;
  int flag = CVode(cvode_mem_, t_out, y_, &t, CV_NORMAL);
  if (flag < 0) { std::fprintf(stderr,"CVode failed, flag=%d\n", flag); exit(EXIT_FAILURE); }
}

