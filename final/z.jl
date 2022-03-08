

using LinearAlgebra, FFTW
import BSON, Statistics, PyPlot, Random

import NMRSpectraSimulator

# for loading something with Interpolations.jl
import OffsetArrays
import Interpolations

#import Clustering

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

load_path = "/home/roy/MEGAsync/inputs/NMR/debug/test_As.bson"
dict = BSON.load(load_path)
As = collect( dict[:As][i] for i = 1:length(dict[:As]) )
