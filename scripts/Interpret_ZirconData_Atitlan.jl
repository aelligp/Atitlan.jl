using JLD2, MagmaThermoKinematics, MAT, CSV, Tables, Statistics, StatsBase
using CairoMakie


function compute_zircons(directory::String)
    dirname = directory # to be updated to Simulation dir e.g., "Atitlan_MTK_1"


    filename = dirname*"/Atitlan_Zircon_analyses.jld2";
    Tracers, Dikes, time_props, Num = JLD2.load(filename, "Tracers", "Dikes", "time_props", "Num")
    println("Loaded data from: $(filename)\n")

    filename = dirname*"/Atitlan_MTK_Final_Setup";
    println("Loading final setup from: $(filename)")

    Atitlan_MTK_Final = load_GMG(filename);

    Time_vec = time_props.Time_vec; # time in seconds
    Melt_Time = time_props.MeltFraction; # time in seconds
    Tav_magma_Time = time_props.Tav_magma; # time in seconds

    # Time_vec,Melt_Time,Tav_magma_Time, Tav_3D_magma_Time, Tracers, Tracers_grid, Tnew, Phases = JLD2.load(filename, "Time_vec", "Melt_Time", "Tav_magma_Time","Tav_3D_magma_Time","Tracers","Tracers_grid", "Tnew_cpu","Phases_float")

    @show dirname

    useTracersOnGrid = true

    Write_R_input = true

    # SecYear = 3600*24*365.25;
    if Write_R_input

        # Prepare the input for the R-script of Gregor Weber
        # if !useTracersOnGrid
            # Advected tracers
            time_vec     = Tracers.time_vec*1e6*SecYear;
            T_vec        = Tracers.T_vec;
            coord        = Tracers.coord;
            Phi          = Tracers.Phi;
            T            = Tracers.T;
        # else

            # # The ones that stayed on the grid points
            # time_vec     = Tracers_grid.time_vec*1e6*SecYear;
            # T_vec        = Tracers_grid.T_vec;
            # coord        = Tracers_grid.coord;
            # Phi          = Tracers_grid.Phi;
            # T            = Tracers_grid.T;
        # end

        # Sample tracers, based on their probability
        # This is done, as the axisymmetric geometry should be expanded to a 3D volume from which melt would be extracted.
        # This implies that the sampling probability scales with 2πR, with
        Prob_Tracers     = [coord[i][1] for i = 1:length(coord)]; Prob_Tracers = Prob_Tracers./maximum(Prob_Tracers)

        if 1==1
            # All Tracers that have final T>700

            ind      = [ T_vec[i][end]>700 for i=1:length(T_vec) ];     # Tracers that are still molten @ the end

            # All Tracers that are still partially molten at the end
            #ind         = findall(Phi .> 0.0)

            time_vec     = time_vec[ind]
            T_vec        = T_vec[ind]
            Prob_Tracers = Prob_Tracers[ind]
            T = T[ind]
        end

        SampleTracers = true
        nSample = 5000
        if SampleTracers
            # Sample the tracers with the probability that they are being extracted:
            id = 1:length(Prob_Tracers)
            # wt = StatsBase.ProbabilityWeights(Prob_Tracers, sum(Prob_Tracers))
            wt = StatsBase.ProbabilityWeights(Prob_Tracers)

            indSample = zeros(Int64,nSample)
            for i=1:nSample
                indSample[i] = sample(id,wt)
            end

            time_vec     = time_vec[indSample]
            T_vec        = T_vec[indSample]
            Prob_Tracers = Prob_Tracers[indSample]
        end

        if 1==1
            # Only consider the most recent time the tracer was molten
            for iT = 1:nSample
                ind             = findall(T_vec[iT] .< 700) #
                if length(ind)>0

                    T_vec[iT]       = T_vec[iT][maximum(ind):end]
                    time_vec[iT]    = time_vec[iT][maximum(ind):end]
                end
            end
        end

        if 1==1
            filename_small = dirname*"/Atitlan_Zircon_analyses_small.jld2";

            # Save it into a different (smaller) file
            jldsave(filename_small; Time_vec, Melt_Time, Tav_magma_Time, time_vec, T_vec,Prob_Tracers)


            #load with:
        # Time_vec,Melt_Time,Tav_magma_Time, time_vec, T_vec, Prob_Tracers = JLD2.load(filename_small, "Time_vec", "Melt_Time", "Tav_magma_Time","time_vec","T_vec","Prob_Tracers")

        end

        # Convert T-t path as text-file input to R-script:
        time_sec, Tt_paths_Temp = compute_zircons_convert_vecs2mat(time_vec, T_vec);

        # if useTracersOnGrid
        #     Tracer_str = "TracersOnGrid";
        # else
            Tracer_str = "Tracers";
        # end


        time_sec    = range(round(time_sec[1]),round(time_sec[end]),length(time_sec));       # avoids round-off errors

        step        = 1;
        Data_TXT                = zeros(length(time_sec), length(1:step:size(Tt_paths_Temp,2))+1);

        Data_TXT[:,1]           = time_sec;

        Data_TXT[:,2:end]       = Tt_paths_Temp[:,1:step:end];

        CSV.write("$(dirname)/$(dirname)_$(Tracer_str).txt", Tables.table(Data_TXT),writeheader=false)

        # Compute average T of points that have T>0
        #  note that the likelihood that something is sampled in 3D is already taken care off above
        T_average   = [ mean(Tt_paths_Temp[i,Tt_paths_Temp[i,:].>0]) for i=1:size(Tt_paths_Temp,1) ]
        time_years  = time_sec/SecYear
        time_Ma     = (time_years[end] .- time_years)/1e6

        # save to CSV file
        ArrayOut = zeros(length(time_Ma),4);
        ArrayOut[:,1] = time_Ma;
        ArrayOut[:,2] = T_average;
        CSV.write("$(dirname)/$(dirname)_$(Tracer_str)_Taverage_julia.csv", Tables.table(ArrayOut),writeheader=false)

        #Plot average T

        fig = Figure(fontsize = 25, size = (800,800))
        ax  = Axis(fig[1, 1],
            xlabel="Age [ka]",
            ylabel="T_average magma [ᵒC]")
        CairoMakie.lines!(ax,time_Ma*1000,T_average, linewidth=0.5)
        CairoMakie.vlines!(ax, 158.0, 0, 100, color=:red, label="W-eruption", linewidth=2)
        CairoMakie.vlines!(ax, 75.0, 0, 100, color=:blue, label="LCY-eruption", linewidth=2)
        CairoMakie.vlines!(ax, 56.0, 0, 100, color=:green, label="I-eruption", linewidth=2)
        axislegend(ax, position=:lt,labelsize=15)
        #limits!(ax,0,1500,600,1000)
        save("$(dirname)/$(dirname)_Taverage_$(Tracer_str).png", fig)


    end

    # use our julia implementation to compute the zircon age distribution
    Compute_with_julia = true
    if Compute_with_julia
    # perform the same zircon age calculations but using julia

    if 1==0
        filename_small = dirname*"/Atitlan_Zircon_analyses_small.jld2";
        # Load 'small' file
        Time_vec, Melt_Time, Tav_magma_Time,time_vec,T_vec,Prob_Tracers = JLD2.load(filename_small, "Time_vec", "Melt_Time", "Tav_magma_Time","time_vec","T_vec","Prob_Tracers")
    end

    # Check if we have valid data before proceeding
    if isempty(time_vec) || isempty(T_vec)
        println("Warning: No valid tracer data found for zircon analysis")
        println("Number of tracers: $(length(time_vec))")
        return  # Exit early if no data
    end

    # Check if any T_vec entries are empty
    valid_tracers = [!isempty(T_vec[i]) for i in 1:length(T_vec)]
    if !any(valid_tracers)
        println("Warning: All tracer temperature vectors are empty")
        return
    end

    # Filter out empty temperature vectors
    time_vec = time_vec[valid_tracers]
    T_vec = T_vec[valid_tracers]

    println("Processing $(length(time_vec)) valid tracers for zircon analysis")

    ZirconData  	=   ZirconAgeData(Tsat=820, Tmin=700, Tsol=700, Tcal_max=800, Tcal_step=1.0, max_x_zr=0.001, zircon_number=100);	 # data as used in the R-script of Gregor

    # time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time, cumPDF = compute_zircons_Ttpath(time_vec/SecYear, T_vec, ZirconData=ZirconData)
    time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average, time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time, cumPDF = compute_zircon_age_PDF(time_vec/SecYear, T_vec, ZirconData=ZirconData, bandwidth = 5.0e4, n_analyses = 500)

        zircon_cumulativePDF = (1.0 .- cumsum(prob))*100;

        # The large numbers in seconds can cause roundoff errors: make sure it is equally spaced again:
        #time_yrs = range(round(time_sec[1]/SecYear),round(time_sec[end]/SecYear),length(time_sec));
        #time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time, zircon_cumulativePDF = compute_zircons_Ttpath(Vector(time_yrs), Tt_paths_Temp, ZirconData=ZirconData)

        ArrayOut = zeros(length(time_sec),5);

        ArrayOut[:,1] = time_sec/SecYear;
        ArrayOut[:,2] = T_av_time;
        ArrayOut[:,3] = number_zircons;
        ArrayOut[:,4] = prob;
        ArrayOut[:,5] = zircon_cumulativePDF;

        # save output to the same directory
        CSV.write("$(dirname)/ZirconPDF_Tav_$(Tracer_str)_$(dirname)_julia.csv", Tables.table(ArrayOut),writeheader=true,header=["time[yrs]","T_average[C]","#zircons","probability","cumulativePDF"])
        time_Ma_vec = (time_years[end] .- time_years)/1e6;

        fig = Figure(fontsize = 25, size = (800,1600))
        ax2  = Axis(fig[1, 1],
        xlabel="Age [ka]",
        ylabel="cumulative PDF ")
        CairoMakie.lines!(ax2,time_Ma_vec*1e3,zircon_cumulativePDF[:], linewidth=2)
        CairoMakie.vlines!(ax2, 158.0, 0, 100, color=:red, label="W-eruption", linewidth=2)
        CairoMakie.vlines!(ax2, 75.0, 0, 100, color=:blue, label="LCY-eruption", linewidth=2)
        CairoMakie.vlines!(ax2, 56.0, 0, 100, color=:green, label="I-eruption", linewidth=2)
        axislegend(ax2, position=:lt,labelsize=15)


        # limits!(ax2,0,(Num.time)/(3600*24*365),0,100)
        ax  = Axis(fig[2, 1],
            xlabel="Age [ka]",
            ylabel="T_average magma [ᵒC]")
        CairoMakie.lines!(ax,time_Ma_vec*1000,T_average, linewidth=0.5, label="from tracers")
        CairoMakie.vlines!(ax, 158.0, 0, 100, color=:red, label="W-eruption", linewidth=2)
        CairoMakie.vlines!(ax, 75.0, 0, 100, color=:blue, label="LCY-eruption", linewidth=2)
        CairoMakie.vlines!(ax, 56.0, 0, 100, color=:green, label="I-eruption", linewidth=2)
        # Time_Vec_ka =  (Time_vec[end] .- Time_vec[:])/SecYear/1e3
        # CairoMakie.lines!(ax,Time_Vec_ka,Tav_3D_magma_Time[:], linewidth=0.5, label="from grid")

        axislegend(ax, position=:lt,labelsize=15)


        # limits!(ax,0,(Num.time)/(3600*24*365),600,1000)

        save("$(dirname)/Tav_Zircon_cumPDF_$(dirname)_$(Tracer_str).png", fig)

        # time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average = zircon_age_PDF(ages_eruptible, number_zircons;  bandwidth=1e4, n_analyses=500, ZirconData=ZirconData)
        #Taken from GeoParams and adapted to our needs
        f = Figure()
        ax = Axis(f[1, 1], xlabel = "Age [Kyr]", ylabel = "Kernel density [ ]", title = "Zircon age probability distribution")
        for i in 1:length(PDF_zircons)
            CairoMakie.lines!((time_Ma[i][end] .- time_Ma[i]) / 1.0e3, PDF_zircons[i], color = "gray66", linewidth = 0.25)
        end
        CairoMakie.lines!(ax, (time_Ma_average[end] .- time_Ma_average) / 1.0e3, PDF_zircon_average, color = "grey0", linewidth = 2.0)
        CairoMakie.xlims!(0.0, 0.4e6 / 1.0e3)
         CairoMakie.vlines!(ax, 158.0, 0, 100, color=:red, label="W-eruption", linewidth=2)
        CairoMakie.vlines!(ax, 75.0, 0, 100, color=:blue, label="LCY-eruption", linewidth=2)
        CairoMakie.vlines!(ax, 56.0, 0, 100, color=:green, label="I-eruption", linewidth=2)
        axislegend(ax, position=:rt, labelsize=15)
        save("$(dirname)/ZirconAgePDF_$(dirname)_$(Tracer_str).png", f)



        # Save data as CSV files
        ArrayOut = zeros(length(time_sec),2);
        ArrayOut[:,1] = time_sec/SecYear;
        ArrayOut[:,2] = T_av_time;
        CSV.write("$(dirname)/$(dirname)_$(Tracer_str)_julia_Taverage.csv", Tables.table(ArrayOut),writeheader=true,header=["time[yrs]","T_average_all_above700[C]"])


        ArrayOut = zeros(length(time_sec),2);
        ArrayOut[:,1] = time_sec/SecYear;
        ArrayOut[:,2] = number_zircons;
        CSV.write("$(dirname)/$(dirname)_$(Tracer_str)_julia_ZirconAges.csv", Tables.table(ArrayOut),writeheader=true,header=["time[yrs]","#zircons"])


    end

end
