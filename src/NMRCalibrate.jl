module NMRCalibrate

import BSON
import NMRData, NMRDataSetup, NMRSpectraSimulator

# Write your package code here.
include("../src/IO/synthetic_setup.jl")

include("../src/DSP/dsp.jl")

export getwraparoundDFTfreqs

end
