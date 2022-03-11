
import NMRDataSetup
import NMRSpectraSimulator
import NMRSpecifyRegions

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

using LinearAlgebra
using FFTW
import PyPlot
import BSON
import JLD

#import Clustering
import Statistics
import Random

include("./helpers/final_helpers.jl")
include("./helpers/solute_calibrate.jl")
include("./helpers/loop_entries1.jl")

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


p_set = Vector{String}(undef, 0)
m_set = Vector{Vector{String}}(undef, 0)
w_set = Vector{Vector{Float64}}(undef, 0)

setupcalibrateserine!(p_set, m_set, w_set)
#setupcalibrateentries1!(p_set, m_set, w_set)
#setupcalibrateentriesNamJan2022!(p_set, m_set, w_set)

projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final"
base_path_JLD = "/home/roy/Documents/repo/NMRData//src/input/molecules"
cs_config_path = "/home/roy/Documents/repo/NMRData/src/input/reduced_cs_config.txt"

function loopscript(p_name_set, m_names_set, w_set,
    max_iters, projects_dir, base_path_JLD, cs_config_path)

    for i = 1:length(p_name_set)
        project_name = p_name_set[i]
        molecule_names = m_names_set[i]
        w = w_set[i]

        calibratesolute(project_name, molecule_names, w;
            max_iters = max_iters
            projects_dir = projects_dir,
            base_path_JLD = base_path_JLD,
            cs_config_path = cs_config_path)

        println()
    end

end


### batch.
max_iters = 50000
loopscript(p_set, m_set, w_set, max_iters, projects_dir, base_path_JLD, cs_config_path)
### end batch.

# #### singular.
# # project_name = "Nam2022_Serine"
# # molecule_names = ["D-(+)-Glucose"; "L-Serine";]
# # w = [1.0; 1.0] #
#
# println("Now on $(project_name)")
# max_iters = 50000
# #max_iters = 5
# include("solute_calibrate.jl")
# println()
# ### end singular.
