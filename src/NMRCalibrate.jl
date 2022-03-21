module NMRCalibrate

#import JLD

using LinearAlgebra
import NLopt
import BSON, Statistics

#import NMRData
import NMRDataSetup # https://github.com/AI4DBiological-Systems/NMRDataSetup.jl
import BoundedLeastSquares # https://github.com/nboyd/BoundedLeastSquares.jl

# dependencies to NMRSpectraSimulator. Need to add explicitly since these packages are not on the Julia public registry.
import GISSMOReader # https://github.com/AI4DBiological-Systems/GISSMOReader.jl

# dependencies to NMRSpecifyRegions. Need to add explicitly since these packages are not on the Julia public registry.
import NMRSpectraSimulator # https://github.com/AI4DBiological-Systems/NMRSpectraSimulator.jl

#
import NMRSpecifyRegions # https://github.com/AI4DBiological-Systems/NMRSpecifyRegions


include("../src/types.jl")
include("../src/utils.jl")

include("../src/IO/synthetic_setup.jl")

include("../src/DSP/dsp.jl")

include("../src/cost/LS.jl")
include("../src/cost/LS_z.jl")

include("../src/cost/updates.jl")
include("../src/cost/multiplier_updates.jl")

include("../src/cost/prep.jl")
include("../src/front_end.jl")

include("../src/cost/nested.jl")

export getwraparoundDFTfreqs

end
