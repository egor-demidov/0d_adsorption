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


// ---------- User data ----------
struct UserData {
  // Parameters to fit:
  sunrealtype k_ads;   // k_a
  sunrealtype k_des;   // k_de
  sunrealtype k_rxn;   // k_r
  sunrealtype S_tot;   // S
  sunrealtype P_tot;   // P_T

  // Known constants:
  sunrealtype k_diff;  // k_d (known)
  sunrealtype a;       // F/V
  sunrealtype X_in;    // inlet
  sunrealtype R;       // radius (constant)
  sunrealtype t_ads_start;
  sunrealtype t_ads_end;
  sunrealtype k_ads_smooth;
};

// Parameter count
static constexpr int NP = 5;
static constexpr int NX = 4;
static constexpr int N  = NX + NP*NX; // 24

// Convenience: get sensitivity block base index
inline int sbase(int j) { return NX + NX*j; }

// ---------- RHS for augmented system ----------
static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  auto* ud = static_cast<UserData*>(user_data);

  const sunrealtype k_a  = ud->k_ads;
  const sunrealtype k_de = ud->k_des;
  const sunrealtype k_r  = ud->k_rxn;
  const sunrealtype S    = ud->S_tot;
  const sunrealtype PT   = ud->P_tot;

  const sunrealtype k_d  = ud->k_diff;
  const sunrealtype a    = ud->a;
  const sunrealtype Xin  = ud->X_in;
  const sunrealtype R    = ud->R;
  const sunrealtype k_a_eff = 0.5 * k_a * (tanh((t - ud->t_ads_start) / ud->k_ads_smooth) - tanh((t - ud->t_ads_end) / ud->k_ads_smooth));

  sunrealtype* Y  = N_VGetArrayPointer(y);
  sunrealtype* dY = N_VGetArrayPointer(ydot);

  // States
  const sunrealtype X   = Y[0];
  const sunrealtype Xgs = Y[1];
  const sunrealtype Xs  = Y[2];
  const sunrealtype P   = Y[3];

  // ---- f(x,theta) ----
  const sunrealtype f1 = a*(Xin - X) - k_d*(X - Xgs);
  const sunrealtype f2 = k_d*(X - Xgs) - (2.0/R)*k_a_eff*Xgs*(S - Xs) + k_de*Xs;
  const sunrealtype f3 = k_a_eff*Xgs*(S - Xs) - (R/2.0)*k_de*Xs - k_r*Xs*(PT - P);
  const sunrealtype f4 = k_r*Xs*(PT - P);

  dY[0] = f1;
  dY[1] = f2;
  dY[2] = f3;
  dY[3] = f4;

  // ---- A = df/dx (4x4) ----
  const sunrealtype A11 = -a - k_d;
  const sunrealtype A12 =  k_d;

  const sunrealtype A21 =  k_d;
  const sunrealtype A22 = -k_d - (2.0/R)*k_a_eff*(S - Xs);
  const sunrealtype A23 =  (2.0/R)*k_a_eff*Xgs + k_de;

  const sunrealtype A32 =  k_a_eff*(S - Xs);
  const sunrealtype A33 = -k_a_eff*Xgs - (R/2.0)*k_de - k_r*(PT - P);
  const sunrealtype A34 =  k_r*Xs;

  const sunrealtype A43 =  k_r*(PT - P);
  const sunrealtype A44 = -k_r*Xs;

  // For each parameter, build b^(j) = df/dtheta_j and compute sdot = A*s + b
  for (int j = 0; j < NP; ++j) {
    const int b = sbase(j);
    const sunrealtype sX   = Y[b+0];
    const sunrealtype sXgs = Y[b+1];
    const sunrealtype sXs  = Y[b+2];
    const sunrealtype sP   = Y[b+3];

    // A*s
    sunrealtype As0 = A11*sX + A12*sXgs;
    sunrealtype As1 = A21*sX + A22*sXgs + A23*sXs;
    sunrealtype As2 = A32*sXgs + A33*sXs + A34*sP;
    sunrealtype As3 = A43*sXs + A44*sP;

    // b^(j)
    sunrealtype b0=0, b1=0, b2=0, b3=0;

    switch (j) {
      case 0: // k_ads
        b1 = -(2.0/R)*Xgs*(S - Xs);
        b2 =  (1.0)*Xgs*(S - Xs);
        break;
      case 1: // k_des
        b1 =  Xs;
        b2 = -(R/2.0)*Xs;
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

  return 0;
}

// ---------- Analytic Jacobian J = d(rhs)/dy (24x24) ----------
// Dense Jacobian signature for CVODE:
// int JacFn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
//           void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  (void)t; (void)fy; (void)tmp1; (void)tmp2; (void)tmp3;
  auto* ud = static_cast<UserData*>(user_data);

  const sunrealtype k_a  = ud->k_ads;
  const sunrealtype k_de = ud->k_des;
  const sunrealtype k_r  = ud->k_rxn;
  const sunrealtype S    = ud->S_tot;
  const sunrealtype PT   = ud->P_tot;

  const sunrealtype k_d  = ud->k_diff;
  const sunrealtype a    = ud->a;
  const sunrealtype R    = ud->R;
  const sunrealtype k_a_eff = 0.5 * k_a * (tanh((t - ud->t_ads_start) / ud->k_ads_smooth) - tanh((t - ud->t_ads_end) / ud->k_ads_smooth));

  sunrealtype* Y = N_VGetArrayPointer(y);

  const sunrealtype X   = Y[0];
  const sunrealtype Xgs = Y[1];
  const sunrealtype Xs  = Y[2];
  const sunrealtype P   = Y[3];

  // ---- A = df/dx ----
  const sunrealtype A11 = -a - k_d;
  const sunrealtype A12 =  k_d;

  const sunrealtype A21 =  k_d;
  const sunrealtype A22 = -k_d - (2.0/R)*k_a_eff*(S - Xs);
  const sunrealtype A23 =  (2.0/R)*k_a_eff*Xgs + k_de;

  const sunrealtype A32 =  k_a_eff*(S - Xs);
  const sunrealtype A33 = -k_a_eff*Xgs - (R/2.0)*k_de - k_r*(PT - P);
  const sunrealtype A34 =  k_r*Xs;

  const sunrealtype A43 =  k_r*(PT - P);
  const sunrealtype A44 = -k_r*Xs;

  // Helper to set dense matrix entry
  auto setJ = [&](int row, int col, sunrealtype val) {
    SM_ELEMENT_D(J, row, col) = val;
  };

  // Zero J first (24x24 is tiny; simplest/robust)
  for (int r = 0; r < N; ++r)
    for (int c = 0; c < N; ++c)
      setJ(r,c,0.0);

  // Top-left block: A
  setJ(0,0, A11); setJ(0,1, A12);
  setJ(1,0, A21); setJ(1,1, A22); setJ(1,2, A23);
  setJ(2,1, A32); setJ(2,2, A33); setJ(2,3, A34);
  setJ(3,2, A43); setJ(3,3, A44);

  // For each sensitivity block j:
  // diagonal block = A, and left block = D^(j) = d/dx (A*s + b^(j))
  for (int j = 0; j < NP; ++j) {
    const int b = sbase(j);

    const sunrealtype sX   = Y[b+0];
    const sunrealtype sXgs = Y[b+1];
    const sunrealtype sXs  = Y[b+2];
    const sunrealtype sP   = Y[b+3];

    // ---- Put diagonal A into block (b..b+3, b..b+3) ----
    setJ(b+0, b+0, A11); setJ(b+0, b+1, A12);

    setJ(b+1, b+0, A21); setJ(b+1, b+1, A22); setJ(b+1, b+2, A23);

    setJ(b+2, b+1, A32); setJ(b+2, b+2, A33); setJ(b+2, b+3, A34);

    setJ(b+3, b+2, A43); setJ(b+3, b+3, A44);

    // ---- Build D^(j) columns via:
    // col_X     = 0
    // col_Xgs   = (dA/dXgs)*s + db/dXgs
    // col_Xs    = (dA/dXs)*s  + db/dXs
    // col_P     = (dA/dP)*s   + db/dP

    // (dA/dXgs)*s = [0, (2/R)k_a*sXs, -k_a*sXs, 0]^T
    const sunrealtype dA_Xgs_0 = 0.0;
    const sunrealtype dA_Xgs_1 = (2.0/R)*k_a_eff*sXs;
    const sunrealtype dA_Xgs_2 = -k_a_eff*sXs;
    const sunrealtype dA_Xgs_3 = 0.0;

    // (dA/dXs)*s = [0, (2/R)k_a*sXgs, -k_a*sXgs, 0]^T
    const sunrealtype dA_Xs_0 = 0.0;
    const sunrealtype dA_Xs_1 = (2.0/R)*k_a_eff*sXgs;
    const sunrealtype dA_Xs_2 = -k_a_eff*sXgs;
    const sunrealtype dA_Xs_3 = 0.0;

    // (dA/dP)*s = [0, 0, k_r*sXs, -k_r*sXs]^T
    const sunrealtype dA_P_0 = 0.0;
    const sunrealtype dA_P_1 = 0.0;
    const sunrealtype dA_P_2 = k_r*sXs;
    const sunrealtype dA_P_3 = -k_r*sXs;

    // db/d(state) depends on j:
    // initialize to 0
    sunrealtype db_Xgs_0=0, db_Xgs_1=0, db_Xgs_2=0, db_Xgs_3=0;
    sunrealtype db_Xs_0 =0, db_Xs_1 =0, db_Xs_2 =0, db_Xs_3 =0;
    sunrealtype db_P_0  =0, db_P_1  =0, db_P_2  =0, db_P_3  =0;

    switch (j) {
      case 0: // k_ads: b = [0, -(2/R)Xgs(S-Xs), Xgs(S-Xs), 0]
        db_Xgs_1 = -(2.0/R)*(S - Xs);
        db_Xgs_2 =  (S - Xs);
        db_Xs_1  =  (2.0/R)*Xgs;
        db_Xs_2  = -Xgs;
        break;
      case 1: // k_des: b = [0, Xs, -(R/2)Xs, 0]
        db_Xs_1  =  1.0;
        db_Xs_2  = -(R/2.0);
        break;
      case 2: // k_rxn: b = [0,0, -Xs(PT-P), +Xs(PT-P)]
        db_Xs_2  = -(PT - P);
        db_Xs_3  =  (PT - P);
        db_P_2   =  Xs;
        db_P_3   = -Xs;
        break;
      case 3: // S_tot: b = [0, -(2/R)k_a Xgs, k_a Xgs, 0]
        db_Xgs_1 = -(2.0/R)*k_a_eff;
        db_Xgs_2 =  k_a_eff;
        break;
      case 4: // P_tot: b = [0,0, -k_r Xs, +k_r Xs]
        db_Xs_2  = -k_r;
        db_Xs_3  =  k_r;
        break;
      default:
        break;
    }

    // Now assemble D^(j) columns (4x4). Column order is [X, Xgs, Xs, P].
    // Column X is zero -> no need to set (already zero).

    // Column Xgs:
    setJ(b+0, 1, dA_Xgs_0 + db_Xgs_0);
    setJ(b+1, 1, dA_Xgs_1 + db_Xgs_1);
    setJ(b+2, 1, dA_Xgs_2 + db_Xgs_2);
    setJ(b+3, 1, dA_Xgs_3 + db_Xgs_3);

    // Column Xs:
    setJ(b+0, 2, dA_Xs_0 + db_Xs_0);
    setJ(b+1, 2, dA_Xs_1 + db_Xs_1);
    setJ(b+2, 2, dA_Xs_2 + db_Xs_2);
    setJ(b+3, 2, dA_Xs_3 + db_Xs_3);

    // Column P:
    setJ(b+0, 3, dA_P_0 + db_P_0);
    setJ(b+1, 3, dA_P_1 + db_P_1);
    setJ(b+2, 3, dA_P_2 + db_P_2);
    setJ(b+3, 3, dA_P_3 + db_P_3);

    // NOTE: these entries fill D^(j) into columns of the *state* block (0..3).
    // We used columns 1,2,3 explicitly; column 0 remains zero. That matches:
    // D^(j)_{:,X}=0, D^(j)_{:,Xgs}, D^(j)_{:,Xs}, D^(j)_{:,P}.
  }

  return 0;
}

// ---------- Minimal CVODE setup example ----------
int main()
{
  // --- user data (example numbers; replace with yours) ---

  const double Di = 40.0;
  const double R = 1.56 / 2.0;
  const double L = 2.0;
  const double V = M_PI * R * R * L;
  const double F = 2.493333333333333 * 760.0 / 1.965;
  const double X_in = 29424160260.695;

  UserData ud{};
  ud.k_ads  = 3.827187367718609e-12;
  ud.k_des  = 0.04253902879091293;
  ud.k_rxn  = 2.0493834170301142e-16;
  ud.S_tot  = 34482079719017.1;
  ud.P_tot  = 88448995701717.02;

  ud.k_diff = 3.66 * Di / R / R;     // known
  ud.a      = F / V;     // F/V
  ud.X_in   = X_in;
  ud.R      = R;

  ud.t_ads_start = 266.3915963445649;
  ud.t_ads_end = 543.5885793335225;
  ud.k_ads_smooth = 2.645834006658198;

  // Initial conditions for x:
  sunrealtype X0=X_in, Xgs0=X_in, Xs0=0, P0=0;

  SUNContext sunctx = nullptr;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx) != 0) {
    std::fprintf(stderr, "SUNContext_Create failed\n");
    return 1;
  }


  // Create N_Vector y0
  N_Vector y = N_VNew_Serial(N, sunctx);
  sunrealtype* Y = N_VGetArrayPointer(y);
  for (int i=0;i<N;++i) Y[i]=0.0;
  Y[0]=X0; Y[1]=Xgs0; Y[2]=Xs0; Y[3]=P0;
  // Sensitivities s(0)=0 already set.

  // Create CVODE memory
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);  // stiff -> BDF
  if (!cvode_mem) { std::fprintf(stderr,"CVodeCreate failed\n"); return 1; }

  sunrealtype t0 = 0.0;
  int flag = CVodeInit(cvode_mem, rhs, t0, y);
  if (flag != CV_SUCCESS) { std::fprintf(stderr,"CVodeInit failed\n"); return 1; }

  CVodeSetMaxNumSteps(cvode_mem, 200000); // or 1000000 if needed

  flag = CVodeSetUserData(cvode_mem, &ud);
  if (flag != CV_SUCCESS) { std::fprintf(stderr,"CVodeSetUserData failed\n"); return 1; }

  // tolerances
  sunrealtype reltol = 1e-8;
  sunrealtype abstol = 1e-10;
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (flag != CV_SUCCESS) { std::fprintf(stderr,"CVodeSStolerances failed\n"); return 1; }

  // Dense linear solver + dense Jacobian
  SUNMatrix A = SUNDenseMatrix(N, N, sunctx);;
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);

  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (flag != CV_SUCCESS) { std::fprintf(stderr,"CVodeSetLinearSolver failed\n"); return 1; }

  flag = CVodeSetJacFn(cvode_mem, jac);
  if (flag != CV_SUCCESS) { std::fprintf(stderr,"CVodeSetJacFn failed\n"); return 1; }

  const long N_t = 800;
  std::vector<double> ts(N_t), Xs(N_t);

  ts[0] = 0.0;
  Xs[0] = X0;
  for (long i = 1; i < N_t; i ++) {
    sunrealtype t;
    const sunrealtype tout =static_cast<double>(i) * 1.0;
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (flag < 0) { std::fprintf(stderr,"CVode failed, flag=%d\n", flag); return 1; }

    ts[i] = tout;

    const sunrealtype *ydata = N_VGetArrayPointer(y);
    Xs[i] = ydata[0];
  }

  std::ofstream ofs("out_data.csv");

  ofs << "t, X\n";
  for (long i = 0; i < N_t; i ++) {
    ofs << ts[i] << ", " << Xs[i] << "\n";
  }

  // Integrate to some time t1
  // sunrealtype t = t0;
  // sunrealtype t1 = 500.0;
  // flag = CVode(cvode_mem, t1, y, &t, CV_NORMAL);
  // if (flag < 0) { std::fprintf(stderr,"CVode failed, flag=%d\n", flag); return 1; }

  // Example: extract Jacobian row for X wrt parameters at t=t1:
  // dX/dtheta = [s^{(0)}_X, s^{(1)}_X, ..., s^{(4)}_X]
  // std::printf("t=%g, X=%g\n", (double)t, (double)Y[0]);
  // for (int j=0;j<NP;++j) {
    // int b = sbase(j);
    // std::printf("  dX/dtheta[%d] = %g\n", j, (double)Y[b+0]);
  // }

  // Cleanup
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);

  return 0;
}
