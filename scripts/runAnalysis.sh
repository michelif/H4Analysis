#!/bin/bash

function createAnalysisTree ()
{
    echo "creating root file for run: $1"
    ./bin/H4Reco cfg/ETHLab_Setup.cfg $1
}

function plotter ()
{
    echo "creating plots for run: $1"
    root -l -b -q "macros/fastPlotter.C(\"$1\")"

}

createAnalysisTree $1
plotter $1