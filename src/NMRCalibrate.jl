module NMRCalibrate


import BSON, Statistics
import NMRData, NMRDataSetup, NMRSpectraSimulator

# https://github.com/nboyd/BoundedLeastSquares.jl
import BoundedLeastSquares

include("../src/types.jl")
include("../src/utils.jl")

include("../src/IO/synthetic_setup.jl")

include("../src/DSP/dsp.jl")

include("../src/cost/LS.jl")
include("../src/cost/updates.jl")


export getwraparoundDFTfreqs

end
