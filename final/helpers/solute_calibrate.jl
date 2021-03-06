
# path to the GISSMO Julia storage folder.
#cs_config_path = "/home/roy/MEGAsync/inputs/NMR/configs/reduced_cs_config.txt"

function calibratesolute(project_name, molecule_names, w;
    max_iters = 50000,
    projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final",
    base_path_JLD = "/home/roy/Documents/repo/NMRData//src/input/molecules",
    cs_config_path = "/home/roy/Documents/repo/NMRData/src/input/reduced_cs_config.txt")


    println("Now on $(project_name)")

    PyPlot.close("all")
    fig_num = 1

    Random.seed!(25)
    PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

    ### user inputs.
    # 0.1% DSS is 0.0046 M = 4.6 mM.
    save_folder_path = joinpath(projects_dir, project_name)
    #@assert 1==2

    #### TODO: load from project name, molecule_names, w from BSON files.
    # find a story to write about.

    # ## reboot. #######
    # project_name = "D-(+)-Glucose-700"
    # molecule_names = ["D-(+)-Glucose"; "DSS"]
    # w = [20.0/4.6; 1.0] # BMRB: DSS is 0.1 % => 4.6 mM

    # project_name = "D-(+)-Glucose-NRC-600"
    # molecule_names = ["D-(+)-Glucose";]
    # w = [1.0; ]

    # project_name = "L-Phenylalanine-700"
    # molecule_names = ["L-Phenylalanine"; "DSS"]
    # w = [20/0.5; 1.0] # BMRB-700 phenylalanine: DSS is 500 micro M.

    # project_name = "L-Glutamine-700"
    # molecule_names = ["L-Glutamine"; "DSS"]
    # w = [20.0/0.5; 1.0] # BMRB: DSS is 500 uM => 0.5 mM

    # project_name = "L-Glutamic acid-700"
    # molecule_names = ["L-Glutamic acid"; "DSS"]
    # w = [20.0/4.6; 1.0] # BMRB: DSS is 0.1% => 4.6 mM

    # project_name = "L-Histidine-700"
    # molecule_names = ["L-Histidine"; "DSS"]
    # w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM
    #
    # project_name = "L-Isoleucine-700"
    # molecule_names = ["L-Isoleucine"; "DSS"]
    # w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM

    # project_name = "L-Serine-700"
    # molecule_names = ["L-Serine"; "DSS"]
    # w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM
    #
    #
    # project_name = "L-Alanine-700"
    # molecule_names = ["L-Alanine"; "DSS"]
    # w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM
    #
    # project_name = "L-Threonine-700"
    # molecule_names = ["L-Threonine"; "DSS"]
    # w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM
    #
    # project_name = "L-Tryptophan-700"
    # molecule_names = ["L-Tryptophan"; "DSS"]
    # w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM
    #
    # project_name = "L-Valine-700"
    # molecule_names = ["L-Valine"; "DSS"]
    # w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM

    w = w ./ norm(w) # since the fit data, y, is normalized.

    # proxy-related.
    tol_coherence = 1e-2
    ??_relative_threshold = 0.05
    #??c_partition_radius = 1e-1 # 17 element partition for glucose spin group 2.
    ??c_partition_radius_candidates = [1e-1; 0.3; 0.5; 0.7; 0.8; 0.9]
    #??c_partition_radius_candidates = [0.9;]
    ??0 = 3.4
    ??cs_max = 0.2 # for proxy.
    ??_??_lb = 0.5
    ??_??_ub = 2.5

    # if a ??c_partition_radius candidate gives max partition sizes less than of
    #   equal to this value, then use that candidate.
    early_exit_part_size = 7

    ### end inputs.



    ## load data.
    load_folder_path = joinpath(projects_dir, project_name)
    load_path = joinpath(load_folder_path, "$(project_name).bson")
    dic = BSON.load(load_path)
    s_t = dic[:s_t]
    fs = dic[:fs]
    SW = dic[:SW]
    ??_0ppm = dic[:??_0ppm]

    # normalize.
    s_t = s_t

    hz2ppmfunc = uu->(uu - ??_0ppm)*SW/fs
    ppm2hzfunc = pp->(??_0ppm + pp*fs/SW)

    offset_Hz = ??_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))

    N = length(s_t)
    DFT_s = fft(s_t)
    U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(N, fs, offset_Hz)

    Z = maximum(abs.(DFT_s))
    y = DFT_s[U_inds] ./ Z


    #
    S_U = DFT_s[U_inds]
    P_y = hz2ppmfunc.(U_y)

    PyPlot.figure(fig_num)
    fig_num += 1

    PyPlot.plot(P_y, real.(S_U), label = "data spectrum")

    PyPlot.legend()
    PyPlot.xlabel("ppm")
    PyPlot.ylabel("real")
    PyPlot.title("data spectra")



    ####### mixture proxy.



    ??cs_max_mixture = collect( ??cs_max for i = 1:length(molecule_names))

    println("Timing: trydiff??cradius()")
    @time mixture_params, ??c_partition_radius = NMRCalibrate.trydiff??cradius(??c_partition_radius_candidates,
        molecule_names, base_path_JLD, ??cs_max_mixture, hz2ppmfunc, ppm2hzfunc,
        fs, SW, ??0, ??_0ppm, early_exit_part_size, ??cs_max, tol_coherence, ??_relative_threshold)
    As = mixture_params



    ??S_ppm = NMRSpectraSimulator.getPsnospininfo(As, hz2ppmfunc)
    ??S_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(??S_ppm))

    println("$(project_name): Partition sizes:")
    display(NMRCalibrate.displaypartitionsizes(As[1]))
    println("??c_partition_radius = ", ??c_partition_radius)
    println()

    #@assert 1==2

    # u_min = ppm2hzfunc(-0.5)
    # u_max = ppm2hzfunc(3.0)

    u_offset = 0.5
    u_min = ppm2hzfunc(??S_ppm_sorted[1] - u_offset)
    u_max = ppm2hzfunc(??S_ppm_sorted[end] + u_offset)

    println("Timing: fitproxies!()")
    @time NMRSpectraSimulator.fitproxies!(As;
    #NMRSpectraSimulator.fitproxiessimple!(As;
    ??_??_lb = ??_??_lb,
    ??_??_ub = ??_??_ub,
    u_min = u_min,
    u_max = u_max,
    ??r = 1.0,
    ????_?? = 0.05)

    #
    # BSON.bson("/home/roy/MEGAsync/inputs/NMR/debug/test_As.bson",
    # As = As)
    # @assert 1==33

    ### cost func.
    combinevectors = NMRSpectraSimulator.combinevectors

    ??sys_cs, y_cost_all, U_cost_all, P_cost_all, exp_info, cost_inds,
    cost_inds_set = NMRCalibrate.prepareoptim(cs_config_path, molecule_names, hz2ppmfunc,
    U_y, y, As; region_min_dist = 0.1)

    # visualize cost regions.
    # sort(P_cost_set[1]) # etc..
    P_cost_set = collect( P_y[cost_inds_set[r]] for r = 1:length(cost_inds_set) )


    ### optim all regions. # 900 secs.
    y_cost = y_cost_all
    U_cost = U_cost_all
    P_cost = P_cost_all

    PyPlot.figure(fig_num)
    fig_num += 1

    PyPlot.plot(P_y, real.(y), label = "data spectrum")
    PyPlot.plot(P_cost, real.(y_cost), "^", label = "positions")

    PyPlot.legend()
    PyPlot.xlabel("ppm")
    PyPlot.ylabel("real")
    PyPlot.title("positions against data spectrum, real part")

    #@assert 1==2

    ??_shifts = NMRSpectraSimulator.combinevectors(??sys_cs)

    ##### set up updates.

    P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
    U = ppm2hzfunc.(P)
    #??S_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.??s) ./ (2*??) ) for A in mixture_params )



    ## parameters that affect qs.
    # A.ss_params.d, A.ss_params.??s_??, A.ss_params.??s_??
    # A.d_singlets, A.??s_singlets, A.??s_singlets, A.??_singlets, A.??0, A.??s_??_singlets
    # purposely perturb ??.

    Es = collect( NMRSpectraSimulator.??CompoundFIDType(As[i]) for i = 1:length(As) )


    ??_lb_default = 0.2
    ??_ub_default = 50.0

    #@assert 1==2

    ## fit model.
    println("Timing: calibrateregions()")
    @time cost_inds_set, p_star_set, ??_BLS_set, d_star_set, ??_star_set, ??_star_set,
    proxies_set = calibrateregions(y, U_y, P_y, cost_inds_set,
    ??_shifts, As, fs, SW, w;
    max_iters = max_iters,
    xtol_rel = 1e-7,
    ftol_rel = 1e-12,
    ??_lb_default = ??_lb_default,
    ??_ub_default = ??_ub_default,
    ??_each_lb = 0.9,
    ??_each_ub = 1.1)

    ### save block.

    save_path = joinpath(save_folder_path, "results_full.bson")
    BSON.bson(save_path,
    p_star_set = p_star_set,
    ??_lb_default = ??_lb_default,
    ??_ub_default = ??_ub_default,
    ??_star_set = ??_BLS_set,
    d_star_set = d_star_set,
    ??_star_set = ??_star_set,
    ??_star_set = ??_star_set,
    w = w,
    # proxy setup-related below.
    ??c_partition_radius = ??c_partition_radius,
    tol_coherence = tol_coherence,
    ??_relative_threshold = ??_relative_threshold,
    ??0 = ??0,
    ??cs_max = ??cs_max,
    ??_??_lb = ??_??_lb,
    ??_??_ub = ??_??_ub,
    # experiement/cost-related below.
    cost_inds_set = cost_inds_set,
    ??sys_cs = ??sys_cs,
    y_cost_all = y_cost_all,
    U_cost_all = U_cost_all,
    P_cost_all = P_cost_all,
    #exp_info = exp_info,
    cost_inds = cost_inds,
    As = As,
    y = y,
    U_y = U_y,
    fs = fs,
    SW = SW,
    ??_0ppm = ??_0ppm)
    ## end save block.


    q_U_set = collect( proxies_set[r].(U) for r = 1:length(proxies_set) )

end
