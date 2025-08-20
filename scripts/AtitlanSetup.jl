# Create Unzen setup
using GeophysicalModelGenerator

# Model setup
println(" --- Generating Setup --- ")

# Topography and project it.

# NOTE: The first time you do this, please set this to true, which will download the topography data from the internet and save it in a file
if !isfile("Data/Topo_cart_Atitlan.jld2")
    using GMT, Statistics
    Topo       =   import_topo(lon = [-91.42, -90.94], lat=[14.46, 14.84], file="@earth_relief_03s.grd")

    proj       =   ProjectionPoint(Lon=-91.18, Lat=14.65)
    Topo_cart  =   convert2CartData(Topo, proj)
    Xt,Yt,Zt   =   xyz_grid(-26:.2:26,-21:.2:21,0)
    Topo_cart  =   project_CartData(CartData(Xt,Yt,Zt,(Zt=Zt,)), Topo, proj)

    save_GMG("Data/Topo_cart_Atitlan", Topo_cart)
end
Topo_cart = load_GMG("Data/Topo_cart_Atitlan")

# Create 3D grid of the region
X,Y,Z       =   xyz_grid(-25:.1:25,-20:.1:20,-25:.1:5)
Data_set3D  =   CartData(X,Y,Z,(Phases=zeros(Int64,size(X)),Temp=zeros(size(X))));       # 3D dataset

# Create 2D cross-section
Nx      =   401;  # resolution in x
Nz      =   201;
Data_2D =   cross_section(Data_set3D, Start=(-25,0), End=(25,0), dims=(Nx, Nz))
Data_2D =   addfield(Data_2D,"FlatCrossSection", flatten_cross_section(Data_2D))
Data_2D =   addfield(Data_2D,"Phases", Int64.(Data_2D.fields.Phases))

# Intersect with topography
Below   =   below_surface(Data_2D, Topo_cart)
Data_2D.fields.Phases[Below] .= 1

# Set Moho
@views Data_2D.fields.Phases[Data_2D.z.val .< -30.0] .= 2

# Set T:
Geotherm = 30
Data_2D.fields.Temp .= -Data_2D.z.val*Geotherm
@views Data_2D.fields.Temp[Data_2D.fields.Temp .<   20.0] .= 20.0
@views Data_2D.fields.Temp[Data_2D.fields.Temp .> 1350.0] .= 1350.0

# Mush reservoir
# add_sphere!(Data_2D.fields.Phases, Data_2D.fields.Temp, Data_2D, cen=(0,0,-10), radius=5, phase=ConstantPhase(2), T=ConstantTemp(800));
add_ellipsoid!(Data_2D.fields.Phases, Data_2D.fields.Temp, Data_2D, cen=(0,0,-6), axes=(8,3,3), phase=ConstantPhase(2), T=ConstantTemp(800));
#add_ellipsoid!(Data_2D.fields.Phases, Data_2D.fields.Temp, Data_2D, cen=(-1,0,-11), axes=(3,3,8), StrikeAngle=225, DipAngle=45, phase=ConstantPhase(5), T=ConstantTemp(1200));
#add_ellipsoid!(Data_2D.fields.Phases, Data_2D.fields.Temp, Data_2D, cen=(-0,0,-23), axes=(8,8,2), StrikeAngle=0, DipAngle=0, phase=ConstantPhase(5), T=ConstantTemp(1200));

write_paraview(Data_2D, "Data/Atitlan_Initial_Setup")
