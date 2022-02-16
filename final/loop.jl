
p_set = Vector{String}(undef, 0)
m_set = Vector{Vector{String}}(undef, 0)
w_set = Vector{Vector{Float64}}(undef, 0)

project_name = "L-Serine-700"
molecule_names = ["L-Serine"; "DSS"]
w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)

project_name = "L-Alanine-700"
molecule_names = ["L-Alanine"; "DSS"]
w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)

project_name = "L-Threonine-700"
molecule_names = ["L-Threonine"; "DSS"]
w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)

project_name = "L-Tryptophan-700"
molecule_names = ["L-Tryptophan"; "DSS"]
w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)

project_name = "L-Valine-700"
molecule_names = ["L-Valine"; "DSS"]
w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)

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

loopscript(p_set, m_set, w_set)
