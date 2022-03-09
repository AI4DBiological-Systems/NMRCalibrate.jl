
p_set = Vector{String}(undef, 0)
m_set = Vector{Vector{String}}(undef, 0)
w_set = Vector{Vector{Float64}}(undef, 0)

#include("./helpers/loop_entries1.jl")
include("./helpers/loop_entries_NamJan2022.jl")

function loopscript(p_name_set, m_names_set, w_set)

    for i = 1:length(p_name_set)
        project_name = p_name_set[i]
        molecule_names = m_names_set[i]
        w = w_set[i]

        println("Now on $(project_name)")
        include("solute_calibrate.jl")
        println()

    end

end


### batch.
loopscript(p_set, m_set, w_set)
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
