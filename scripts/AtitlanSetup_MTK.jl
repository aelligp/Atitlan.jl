# Atitlan setup
using JLD2, Random, CairoMakie
const USE_GPU=false;
if USE_GPU
    using CUDA      # needs to be loaded before loading Parallkel=
end
using ParallelStencil, ParallelStencil.FiniteDifferences2D

using MagmaThermoKinematics
@static if USE_GPU
    environment!(:gpu, Float64, 2)      # initialize parallel stencil in 2D
    CUDA.device!(0)                     # select the GPU you use (starts @ zero)
    @init_parallel_stencil(CUDA, Float64, 2)
else
    environment!(:cpu, Float64, 2)      # initialize parallel stencil in 2D
    @init_parallel_stencil(Threads, Float64, 2)
end
using MagmaThermoKinematics.Diffusion2D # to load AFTER calling environment!()

# Import a few routines, so we can overwrite them below
using MagmaThermoKinematics.MTK_GMG     # Allow overwriting user routines
using StructArrays
include("AtitlanSetup.jl")                # Create Model setup

println(" --- Performing MTK models --- ")

# Overwrite some of the default functions
@static if USE_GPU
    function MTK_GMG.MTK_print_output(Grid::GridData, Num::NumericalParameters, Arrays::NamedTuple, Mat_tup::Tuple, Dikes::DikeParameters)
        println("$(Num.it), Time=$(round(Num.time/Num.SecYear)) yrs; max(T) = $(round(maximum(Arrays.Tnew)))")
        return nothing
    end
else
    function MTK_GMG.MTK_visualize_output(Grid::GridData, Num::NumericalParameters, Arrays::NamedTuple, Mat_tup::Tuple, Dikes::DikeParameters)
        if mod(Num.it,Num.CreateFig_steps)==0
            x_1d        =   Grid.coord1D[1]/1e3;
            z_1d        =   Grid.coord1D[2]/1e3;
            temp_data   =   Array(Arrays.Tnew)'
            ϕ_data      =   Array(Arrays.ϕ)'
            phase_data  =   Float64.(Array(Arrays.Phases))'

            # remove topo on plots
            ind             = findall(phase_data .== 0)
            phase_data[ind] .= NaN
            temp_data[ind]  .= NaN

            t = Num.time/SecYear/1e3;

            fig = Figure(fontsize = 25, size = (800,1600))
            ax = Axis(fig[1,1], aspect =  DataAspect(), title="Temperature t=$(round(t)) kyrs", xlabel="x [km]", ylabel="z [km]", limits=(nothing, (minimum(z_1d), 3)))
            ax2 = Axis(fig[1,2], aspect =  DataAspect(), title="Melt fraction", xlabel="x [km]", ylabel="z [km]", limits=(nothing, (minimum(z_1d), 3)))

            heatmap!(ax,x_1d, z_1d, temp_data', colormap=:viridis)
            heatmap!(ax2,x_1d, z_1d, ϕ_data',    colormap=:lipari)

            display(fig)
        end
        return nothing
    end
end

function MTK_GMG.MTK_print_output(Grid::GridData, Num::NumericalParameters, Arrays::NamedTuple, Mat_tup::Tuple, Dikes::DikeParameters)
    println("$(Num.it), Time=$(round(Num.time/Num.SecYear/1e3, digits=3)) kyrs; max(T) = $(round(maximum(Arrays.Tnew))); mean(T) = $(round(mean(Arrays.Tnew[Arrays.Phases .> 1])))")
    println("Dike center: ", Dikes.Center)
    return nothing
end

"""
    MTK_update_Arrays!(Arrays::NamedTuple, Grid::GridData, Dikes::DikeParameters, Num::NumericalParameters)

Update arrays and structs of the simulation (in case you want to change them during a simulation)
You can use this, for example, to change the size and location of an intruded dike
"""
function MTK_GMG.MTK_update_ArraysStructs!(Arrays::NamedTuple, Grid::GridData, Dikes::DikeParameters, Num::NumericalParameters, Mat_tup::Tuple)

    if Num.AddRandomSills && mod(Num.it,Num.RandomSills_timestep)==0
        # This randomly changes the location and orientation of the sills
        if Num.dim==2
            Loc = [Dikes.W_ran; Dikes.H_ran]
        else
            Loc = [Dikes.W_ran; Dikes.L_ran; Dikes.H_ran]
        end

        # Randomly change location of center of dike/sill
        cen       = (0.0, -6000) .+ rand(-0.5:1e-3:0.5, Num.dim).*Loc;

        Dip       = rand(-Dikes.Dip_ran/2.0    :   0.1:   Dikes.Dip_ran/2.0)
        Strike    = rand(-Dikes.Strike_ran/2.0 :   0.1:   Dikes.Strike_ran/2.0)

        if cen[end]<Dikes.SillsAbove;
            Dip = Dip   + 90.0                                          # Orientation: near-vertical @ depth
        end

        Dikes.Center = cen;
        Dikes.Angle  = [Dip, Strike];
    end

    # change dike injecion interval depending on time
    if (Num.time/SecYear>120e3  &&  Num.time/SecYear<140e3)
        Dikes.InjectionInterval = 1000*SecYear;
    elseif (Num.time/SecYear>140e3  &&  Num.time/SecYear<225e3)
        Dikes.InjectionInterval = 1500*SecYear;
    else
        Dikes.InjectionInterval = 2500*SecYear;
    end
    println("Dike injection interval = $(Dikes.InjectionInterval/SecYear) years ")

    return nothing
end


"""
    Tracers = MTK_inject_dikes(Grid, Num, Arrays, Mat_tup, Dikes, Tracers, Tnew_cpu)

Function that injects dikes once in a while
"""
function MTK_GMG.MTK_inject_dikes(Grid::GridData, Num::NumericalParameters, Arrays::NamedTuple, Mat_tup::Tuple, Dikes::DikeParameters, Tracers::StructVector, Tnew_cpu)

    if floor(Num.time/Dikes.InjectionInterval)> Dikes.dike_inj || mean(Arrays.Tnew[Arrays.Phases .> 1]) .< 700.0
        Dikes.dike_inj      =   floor(Num.time/Dikes.InjectionInterval)                 # Keeps track on what was injected already
        if Num.dim==2
            T_bottom  =   Tnew_cpu[:,1]
        else
            T_bottom  =   Tnew_cpu[:,:,1]
        end
        dike      =   Dike(W=Dikes.W_in, H=Dikes.H_in, Type=Dikes.Type, T=Dikes.T_in_Celsius, Center=Dikes.Center[:],  Angle=Dikes.Angle, Phase=Dikes.DikePhase);               # "Reference" dike with given thickness,radius and T
        Tnew_cpu .=   Array(Arrays.T)

        Tracers, Tnew_cpu,Vol,Dikes.dike_poly, VEL  =   InjectDike(Tracers, Tnew_cpu, Grid.coord1D, dike, Dikes.nTr_dike, dike_poly=Dikes.dike_poly);     # Add dike, move hostrocks

        if Num.flux_bottom_BC==false
            # Keep bottom T constant (advection modifies this)
            if Num.dim==2
                Tnew_cpu[:,1]     .=  T_bottom
            else
                Tnew_cpu[:,:,1]   .=  T_bottom
            end
        end

        Arrays.T           .=   Data.Array(Tnew_cpu)
        Dikes.InjectVol    +=   Vol                                                     # Keep track of injected volume
        Qrate               =   Dikes.InjectVol/Num.time
        Dikes.Qrate_km3_yr  =   Qrate*SecYear/km³
        Qrate_km3_yr_km2    =   Dikes.Qrate_km3_yr/(pi*(Dikes.W_in/2/1e3)^2)
        println("  Added new dike; time=$(Num.time/kyr) kyrs, total injected magma volume = $(Dikes.InjectVol/km³) km³; rate Q= $(Dikes.Qrate_km3_yr) km³yr⁻¹")

        if Num.advect_polygon==true && isempty(Dikes.dike_poly)
            Dikes.dike_poly   =   CreateDikePolygon(dike);            # create dike for the 1th time
        end

        if length(Mat_tup)>1
           PhasesFromTracers!(Array(Arrays.Phases), Grid, Tracers, BackgroundPhase=Dikes.BackgroundPhase, InterpolationMethod="Constant");    # update phases from grid

           # Ensure that we keep the initial phase of the area (host rocks are not deformable)
           if Num.keep_init_RockPhases==true
                Phases      = Array(Arrays.Phases)          # move to CPU
                Phases_init = Array(Arrays.Phases_init)
                for i in eachindex(Phases)
                    if Phases[i] != Dikes.DikePhase
                        Phases[i] = Phases_init[i]
                    end
                end
                Arrays.Phases .= Data.Array(Phases)          # move back to GPU
           end
        end

    end

    return Tracers
end


"""
    MTK_update_TimeDepProps!(time_props::TimeDependentProperties, Grid::GridData, Num::NumericalParameters, Arrays::NamedTuple, Mat_tup::Tuple, Dikes::DikeParameters)

Update time-dependent properties during a simulation
"""
function MTK_GMG.MTK_update_TimeDepProps!(time_props::TimeDependentProperties, Grid::GridData, Num::NumericalParameters, Arrays::NamedTuple, Mat_tup::Tuple, Dikes::DikeParameters)
    push!(time_props.Time_vec,      Num.time);   # time
    push!(time_props.MeltFraction,  sum( Arrays.ϕ)/(Num.Nx*Num.Nz));    # melt fraction

    ind = findall(Arrays.T.>700);
    if ~isempty(ind)
        Tav_magma_Time = sum(Arrays.T[ind])/length(ind)     # average T of part with magma
    else
        Tav_magma_Time = NaN;
    end
    push!(time_props.Tav_magma, Tav_magma_Time);        # average magma T
    push!(time_props.Tmax,      maximum(Arrays.T));     # maximum magma T
    push!(time_props.Qrate,       Dikes.Qrate_km3_yr);  # magma flux
    return nothing
end

# Define a new structure with time-dependent properties
@with_kw mutable struct TimeDepProps1 <: TimeDependentProperties
    Time_vec::Vector{Float64}       = [];           # Center of dike
    MeltFraction::Vector{Float64}   = [];           # Melt fraction over time
    Tav_magma::Vector{Float64}      = [];           # Average magma
    Tmax::Vector{Float64}           = [];           # Max magma temperature
    Tmax_1::Vector{Float64}         = [];           # Another magma temperature vector
    Qrate::Vector{Float64}          = [];           # Magma flux
end

# Define numerical parameters
Num         = NumParam( SimName             =   "Atitlan_August2025",
                        dim                 =   2,
                        Nx                  =   200,
                        Nz                  =   200,
                        maxTime_Myrs        =   0.3,
                        SaveOutput_steps    =   25,
                        CreateFig_steps     =   5,
                        USE_GPU             =   USE_GPU,
                        ω                   =   0.5,
                        AddRandomSills      =   true,
                        RandomSills_timestep=   5);

# dike parameters
Dike_params = DikeParam(Type                    =   "ElasticDike",
                        InjectionInterval_year  =   500,       # flux= 14.9e-6 km3/km2/yr
                        Center                  =   [0.0, -6000],        # Center of dike
                        T_in_Celsius            =   1000.0,      # Temperature of dike
                        W_in                    =   5e3,
                        H_in                    =   250,
                        H_ran                   =   1000,       # height of random injection area
                        W_ran                   =   10000,       # width of random injection area
                        nTr_dike                =   2000,
                        Dip_ran                 =   45,         # angle aroun d which we randomly change the dip
                        DikePhase               =   3,          # phase of dike
                        SillsAbove              =   -10e3       # below this we have dikes; above sills
                )

# Define parameters for the different phases
MatParam     = (SetMaterialParams(Name="Air", Phase=0,
                                Density      = ConstantDensity(ρ=2700kg/m^3),
                                LatentHeat   = ConstantLatentHeat(Q_L=0.0J/kg),
                                Conductivity = ConstantConductivity(k=3Watt/K/m),          # in case we use constant k
                                HeatCapacity = ConstantHeatCapacity(Cp=1000J/kg/K),
                                Melting      = MeltingParam_Caricchi()),
                SetMaterialParams(Name="Crust", Phase=1,
                                Density      = ConstantDensity(ρ=2700kg/m^3),
                                LatentHeat   = ConstantLatentHeat(Q_L=3.13e5J/kg),
                                Conductivity = T_Conductivity_Whittington_parameterised(),   # T-dependent k
                                HeatCapacity = ConstantHeatCapacity(Cp=1000J/kg/K),
                                Melting      = MeltingParam_Caricchi()),
                                # Melting      = MeltingParam_Smooth3rdOrder()),
                SetMaterialParams(Name="Mush", Phase=2,
                                Density      = ConstantDensity(ρ=2700kg/m^3),
                                LatentHeat   = ConstantLatentHeat(Q_L=3.13e5J/kg),
                                Conductivity = T_Conductivity_Whittington_parameterised(),   # T-dependent k
                                HeatCapacity = ConstantHeatCapacity(Cp=1000J/kg/K),
                                Melting      = MeltingParam_Caricchi()),
                                # Melting      = MeltingParam_Smooth3rdOrder()),
                SetMaterialParams(Name="Dikes", Phase=3,
                                Density      = ConstantDensity(ρ=2500kg/m^3),
                                LatentHeat   = ConstantLatentHeat(Q_L=3.13e5J/kg),
                                Conductivity = T_Conductivity_Whittington_parameterised(),   # T-dependent k
                                HeatCapacity = ConstantHeatCapacity(Cp=1000J/kg/K),
                                Melting      = MeltingParam_Caricchi()),
                                # Melting      = MeltingParam_Smooth3rdOrder())
                )

# Call the main code with the specified material parameters

# Call the main code with the specified material parameters
flux  = ((4/3*pi*(Dike_params.W_in/2)*(Dike_params.W_in/2)*(Dike_params.H_in/2))./1e9) / Dike_params.InjectionInterval_year;  # flux in km3/yr
printstyled("You have defined a flux of $(round(flux, digits=7)) km³/yr for the dikes.\n", color=:default, bold=true)
#The call below runs the main code
Grid, Arrays, Tracers, Dikes, time_props = MTK_GeoParams_2D(MatParam, Num, Dike_params, CartData_input=Data_2D, time_props=TimeDepProps1()); # start the main code

println(" --- MTK models finished --- ")
println(" --- Saving results dont quite yet --- ")
# Save final data to disk
Data_set2D_out = Data_2D;
# Data_set2D_out = MTK_GMG.add_data_CartData(Data_set2D_out, "Temperature[C]",  Array(Arrays.Tnew ));   # in MPa
Data_set2D_out = MTK_GMG.add_data_CartData(Data_set2D_out, "Temp",         Array(Arrays.Tnew));
Data_set2D_out = MTK_GMG.add_data_CartData(Data_set2D_out, "Phases",       Array( Float64.(Arrays.Phases)));
Data_set2D_out = MTK_GMG.add_data_CartData(Data_set2D_out, "MeltFraction", Array(Arrays.ϕ));

save_GMG("$(Num.SimName)/Atitlan_MTK_Final_Setup", Data_set2D_out)

filename = "$(Num.SimName)/Atitlan_Zircon_analyses.jld2"
jldopen(filename, "w"; iotype=IOStream) do file
    file["Tracers"] = Tracers
    file["Dikes"] = Dikes
    file["time_props"] = time_props
    file["Num"] = Num
end
println("Saved data to: $(filename)\n")


println("Now running Zircon interpretation script")
include("Interpret_ZirconData_Atitlan.jl")  # Interpret the zircon data
compute_zircons(Num.SimName)  # Call the function to compute zircons
