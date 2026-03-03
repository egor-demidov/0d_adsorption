$buildDir = ".\build-release\"
$installDir = ".\dist\bin"
$sourceDir = ".\"

robocopy "$buildDir" "$installDir" *.dll
robocopy "$buildDir" "$installDir" *.exe
robocopy "$sourceDir" "$installDir" preprocess.py
robocopy "$sourceDir" "$installDir" plot_fitted_curve.py
robocopy "$sourceDir" "$installDir" plot_run_only.py
