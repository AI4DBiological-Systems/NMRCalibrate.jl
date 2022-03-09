
project_name = "L-Glutamine-700"
molecule_names = ["L-Glutamine"; "DSS"]
w = [20.0/0.5; 1.0] # BMRB: DSS is 500 uM => 0.5 mM
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)

project_name = "L-Serine-700"
molecule_names = ["L-Serine"; "DSS"]
w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)

project_name = "D-(+)-Glucose-NRC-600"
molecule_names = ["D-(+)-Glucose";]
w = [1.0; ]
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)

project_name = "Nam2022_Serine"
molecule_names = ["D-(+)-Glucose"; "L-Serine";]
w = [1.0; 1.0] #
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)

project_name = "L-Leucine-500-2mM"
molecule_names = ["L-Leucine";]
w = [1.0]
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)

project_name = "L-Isoleucine-700"
molecule_names = ["L-Isoleucine"; "DSS"]
w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM
push!(p_set, project_name)
push!(m_set, molecule_names)
push!(w_set, w)
