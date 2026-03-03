This program fits experimental uptake curves to the chained CSTR model.
Follow the instructions below to process an example uptake curve.
A Windows PowerShell session is needed to follow this tutorial.


    STEP 0. ENVIRONMENT SET UP
    A Python environment needs to be set up to enable execution of pre/postprocessing scripts. A virtual environment
    is recommended:

        python3 -m venv venv    # Create a virtual environment called 'venv' (only done ONCE)
        .\venv\Scripts\activate # Activate the created virtual environment

    Once a virtual environment has been created and activated, dependencies need to be installed into the environment.
    That can be done with:

         pip install -r requirements.txt    # Install the dependencies (only done ONCE)

     Finally, the folder with scripts and executables needs to be added to Path:

        $env:PATH += ";$((Get-Location).Path)\bin"


    STEP 1. PREPROCESSING
    At this point an experimental uptake curve, as well as some parameters, are read from an Excel workbook by a
    Python script. The script performs drift correction, detects times at which adsorption and desorption occur,
    and generates an input file compatible with the CSTR model in .json format. To preprocess the example, run:

        cd .\example\
        preprocess .\uptake_curve.xlsx --worksheet "NaCl-2" --output nacl-2.json

    A plot window will appear. Following the uptake curve from left to right, the user needs to click on the uptake
    curve at points where is known to be zero. The curve will be shifted up/down accordingly to correct for signal
    drift. Once drift correction is completed, the plot window needs to be closed. An output file in .json format will
    be automatically generated then.


    STEP 2. FITTING
    Inspect the file "nacl-2.json". Some of the values stored in that file were read from "uptake_curve.xlsx", others -
    were populated with default values. Namely, initial guesses for kinetic parameters, "k_ads_smooth" (peak smoothing)
    parameter, and "N_reactors" (number of chained CSTR reactors) are populated with default values. Defaults can be
    changed by editing "preprocess.py" script. Otherwise, these values might need to be tweaked before the fitter run.

    To attempt to fit preprocessed experimental data to the CSTR model, run:

        0d_adsorption_fit_chained.exe .\nacl-2.json

    If the run was successful, "Termination: CONVERGENCE" will be printed on the last line of standard output. If the
    run failed, initial guesses may need to be adjusted.


    STEP 3. PLOTTING THE RESULTS
    After a successful fitting run, a file named "fitted.json" will appear in the working directory. It will contain the
    values of fitted parameters, their estimated standard error, as well as uptake curves at the initial guess (X0),
    at the solution (X), and sensitivities. The solution can be quickly visualized and evaluated for correctness by
    running:

        plot_fitted_curve --input nacl-2.json --solution fitted.json

    If the "k_ads_smooth" parameter is not appropriate, peaks on the simulated uptake curve will be too sharp or too
    smooth and the fit will be poor. If that happens, "k_ads_smooth" needs to be adjusted and the curve - refitted until
    the sharpness of experimental peaks is matched.


    STEP 4. PREDICTIVE RUNS
    With kinetic parameters known for a certain reactant - substrate pair, uptake curves at different conditions can be
    predicted. A different executable that only runs the CSTR model once without fitting data needs to be used for
    predictive runs. Input for the predictive model is similar in format to input for the fitting model, but some
    parameters differ. For example, the predictive model needs "t_tot" - the total duration of the run. In the case
    of the fitter, duration is determined automatically based on provided experimental data and "t_tot" does not need
    to be supplied separately. To run the provided example:

        0d_adsorption_chained.exe .\run_only_input.json

    Computed uptake curve is stored in a file named "run.json". To display the generated curve, run:

        plot_run_only .\run.json


