# This script create default configuration files that contain the tuning parameters
# for some of the algorithms this package. We assign a default value for each
# compound name key entry in the user-specified
# "compound to parameters filename" JSON file.

# One should manually edit the generated config files according to
# their preferred tuning parameters.

import JSON
import JSON3

name_map_filename = "select_compounds"
load_path = joinpath("/home/roy/Documents/repo/NMRData/input/compound_mapping", "$(name_map_filename).json")
dict_compound_to_filename = JSON.parsefile(load_path)
compound_names = collect( key for (key,value) in dict_compound_to_filename )


##### config file for SH simulation.
function createdefaultcalibrationconfig(compound_names::Vector{String};
    max_cs_shift::Float64 = 0.05, # in ppm.
    α_relative_threshold::Float64 = 0.05,
    Δc_partition_radius::Float64 = 0.3)

    db_dict = Dict()
    for name in compound_names

        db_dict[name] = Dict("maximum chemical shift" => max_cs_shift,
            "relative amplitude threshold" => α_relative_threshold,
            "maximum Δc deviation" => Δc_partition_radius)
    end

    return db_dict
end

save_folder = "/home/roy/Documents/repo/NMRData/input/SH_configs"
save_name = "$(name_map_filename)_SH_configs_default.json"
save_path = joinpath(save_folder, save_name)

# write the default SH configurations to json.
dict_out = createdefaultSHconfig(compound_names)

stringdata = JSON.json(dict_out)

open(save_path, "w") do f
    JSON3.pretty(f, stringdata)
    println(f)
end
