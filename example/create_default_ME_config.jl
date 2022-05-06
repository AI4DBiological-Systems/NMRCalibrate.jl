
import NMRSpectraSimulator

import Unicode
import JSON
import JSON3
import Destruct

output_folder = "/home/roy/Documents/repo/NMRData/input/fit_configs"
mkpath(output_folder)
save_path = joinpath(output_folder, "calibrate_700MHz_type1_select_compounds_default.json")

# get mapping from molecule names to their spin system info json files.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

# get list of molecule names.
molecule_names = collect( key for (key, value) in dict_compound_to_filename )

# set up machine simulate and get
# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665

println("Timing: setupmixtureproxies()")
@time mixture_params = NMRSpectraSimulator.setupmixtureSH(molecule_names,
    H_params_path, dict_compound_to_filename, fs, SW,
    ν_0ppm)
As = mixture_params

println("Timing: mag equivalence")
@time MEs = NMRSpectraSimulator.getmageqinfomixture(molecule_names,
    H_params_path,
    dict_compound_to_filename;
    unique_cs_atol = 1e-6)

# corresponding to NMRSpectraSimulator.SpinSysFIDType1{T}
function populatedefaultvaluescompound1(N_spin_systems::Int, N_singlets::Int, fill_val::T) where T
    #
    systems = collect( fill_val for i = 1:N_spin_systems)
    singlets = collect( fill_val for i = 1:N_singlets)

    return [systems; singlets]
end


"""
Creates ordering of length N from equivalent indices.
`eq_inds` can only contain unique values from 1:N.
When used for a spin system, N is the number of nuclei in the system.

Returns the indices for the N nuclei, and the number of degrees of freedom.
"""
function createorderingfromeqinds(eq_inds::Vector{Vector{Int}}, N::Int)

    j = 0 # degrees of freedom counter.
    out = zeros(Int, N)

    if isempty(eq_inds)
        return collect(1:N), N
    end

    # check if `eq_inds` only contains unique values from 1:N.
    if all( all( (eq_inds[k][l] in 1:N) for l = 1:length(eq_inds[k])) for k = 1:length(eq_inds) )
        for k = 1:length(eq_inds)
            j += 1

            out[eq_inds[k]] .= j
        end

        # fill the rest.
        for n = 1:N
            if out[n] == 0
                # out[n] is unassigned.
                j += 1
                out[n] = j
            end
        end

        return out, j
    end

    println("Invalid ME, using default.")
    return collect(1:N), N
end

"""
based on the provided default values, and β ordering based on magnetic equivalence of nuclei.
"""
function createcalibrateconfigdict(As, MEs, molecule_names;
    max_cs_shift_default = 0.05,
    λ_lb_default = 0.7,
    λ_ub_default = 1.5)

    @assert length(As) == length(molecule_names) == length(MEs)

    dict_out = Dict()
    for n = 1:length(As)

        ## parse.
        A = As[n]
        name = molecule_names[n]
        ME = MEs[n]

        ## prepare.
        key_out = "$(name)"

        # spin systems and singlets.
        max_cs_shift = populatedefaultvaluescompound1(length(A.N_spins_sys), length(A.αs_singlets), max_cs_shift_default)
        λ_lb = populatedefaultvaluescompound1(length(A.N_spins_sys), length(A.αs_singlets), λ_lb_default)
        λ_ub = populatedefaultvaluescompound1(length(A.N_spins_sys), length(A.αs_singlets), λ_ub_default)

        # β for spin systems only. Local (to each spin system) indices.
        tmp = collect( createorderingfromeqinds(ME[i], A.N_spins_sys[i]) for i = 1:length(A.N_spins_sys) )
        κs_β_ordering, κs_β_DOF = Destruct.destruct(tmp)

        ## store.
        dict_out[key_out] = Dict( "λ_lb" => λ_lb,
        "λ_ub" => λ_ub,
        "maximum chemical shift" => max_cs_shift,
        "κs_β ordering" => κs_β_ordering,
        "κs_β degrees of freedom" => κs_β_DOF)
    end

    return dict_out
end


dict_out = createcalibrateconfigdict(As, MEs, molecule_names;
    max_cs_shift_default = 0.05,
    λ_lb_default = 0.7,
    λ_ub_default = 1.5)

stringdata = JSON.json(dict_out)

open(save_path, "w") do f
    JSON3.pretty(f, stringdata)
    println(f)
end
