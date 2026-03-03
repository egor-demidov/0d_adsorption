#!/bin/bash

BUILD_DIR="./build-release/"
INSTALL_DIR="./dist/bin"
SOURCE_DIR="."

rsync -u "$BUILD_DIR/0d_adsorption_chained" "$INSTALL_DIR"
rsync -u "$BUILD_DIR/0d_adsorption_fit_chained" "$INSTALL_DIR"
rsync -u "$SOURCE_DIR/preprocess.py" "$INSTALL_DIR"
rsync -u "$SOURCE_DIR/plot_fitted_curve.py" "$INSTALL_DIR"
rsync -u "$SOURCE_DIR/plot_run_only.py" "$INSTALL_DIR"
