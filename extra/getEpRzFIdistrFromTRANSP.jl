#########################################  getEpRzFIdistrFromTRANSP.jl #########################################
# This Julia script will read the TRANSP fast-ion (FI) .cdf-file, use read_nubeam() and output a 4D FI-distribution 
# in particle-space coordinates (E,p,R,z), together with the E,p,R and z vectors.
# It will save this data in .jld2 format.
#
# For this to work, you need to manually specify the minimum and maximum values for major radius and vertical
# positions R and z, respectively. You also need to specify how many grid points you want in R and z. These 
# specifications are set via the input arguments below this script description.

### Inputs
# filepath_distr - The path to the TRANSP fast-ion distribution in .cdf file format that you want to use - String
# folderpath_OWCF - The path to the OWCF-folder. Remember to finish with '/' since it's a folder and not a file. - String
# folderpath_out - The path to the folder in which you want your output. Remember to finish with '/'. - String
# nR - The number of major radius position grid points that you would like your (E,p,R,z) distribution to have - Int64
# nz - The number of vertical position grid points that you would like your (E,p,R,z) distribution to have - Int64
# Rmin - The minimum value (lower boundary) for your major radius position grid points. In centimeters. - Float64
# Rmax - The maximum value (lower boundary) for your major radius position grid points. In centimeters. - Float64
# signOfJdotB - Should be set to 1 if sign(dot(J,B))>0. Else, set to -1.
# zmin - The minimum value (lower boundary) for your vertical position grid points. In centimeters. - Float64
# zmax - The maximum value (lower boundary) for your vertical position grid points. In centimeters. - Float64

### Outputs
# - 

### Saved file
## [TRANSP ID]_fi_[step].jld2 with the following keys
# species - The fast-ion species of the fast-ion distribution. Saved as a String
# timepoint - The timepoint of the tokamak shot of the TRANSP simulation. Saved as a Float64
# ntot - The total number of fast ions. Saved as a Float64
# f - The (E,p,R,z) fast-ion distribution. Saved as a 4D matrix
# f_Rz - The f(R,z). The energy and pitch dependence has been integrated out. Saved as a 2D array
# energy - The energy grid points saved as a 1D array
# pitch - The pitch grid points saved as a 1D array
# R - The major radius position grid points saved as a 1D array. cm
# z - The vertical position grid points saved as 1D array. cm
# R_sep - The R coordinates of the separatrix saved as a 1D array. cm
# z_sep - The z coordinates of the separatrix saved as a 1D array. cm
# signOfJdotB - The signOfJdotB input saved as a scalar
# orig_TRANSP_filepath - filepath_distr saved as a string

# Script written by Henrik Järleblad, last updated 2022-10-06.
# The rz_grid() and read_nubeam() functions have been re-written in Julia; the original functions were written and developed in Python
# as part of the FIDASIM code framework, by Luke Stagner.
################################################################################################################

############################################################# Inputs
filepath_distr = "C:/Users/henrikj/Documents/codes/OWCF/TRANSP/JET/99965/K73/99965K73_fi_2.cdf" # for example '/Users/anna/TRANSP/JET/99500/V05/99500V05_fi_1.cdf'
folderpath_OWCF = "G:/My Drive/DTU/codes/OWCF/"
folderpath_out = "C:/Users/henrikj/Documents/codes/OWCF/TRANSP/JET/99965/K73/"
nR = 64
nz = 65
Rmin = 1.83591998*100 # cm. Default value suitable for JET
Rmax = 3.89143991*100 # cm. Default value suitable for JET
zmin = -1.74594998*100 # cm. Default value suitable for JET
zmax = 1.98403001*100 # cm. Default value suitable for JET
e_range = nothing # Specified as (a,b). If all TRANSP energies are desired, leave as 'nothing'
p_range = nothing # Specified as (a,b). If all TRANSP pitch values are desired, leave as 'nothing'
species=1 # '1' means the first species of the TRANSP file (usually the main fast-ion species)
btipsign=1 # '1' means that sign(J ⋅ B) > 0. '-1' means that sign(J ⋅ B) < 0
verbose=true # If true, the script will speak a lot!
plotting=true # If true, the script will plot results along the way. For safety checks.
vverbose=false # If true, the script will speak a lot more! WARNING! You don't want this... It's just too much.
save_plots = true

############################################################# Activate OWCF environment and load packages
using Pkg
cd(folderpath_OWCF)
Pkg.activate(".")

using JLD2
using Plots
using NetCDF
using Statistics
using LinearAlgebra
using Interpolations
using VoronoiDelaunay
using NearestNeighbors
############################################################# Define the grid Struct, rz_grid() and read_nubeam() functions
struct grid
    r2d::Matrix{Float64}
    z2d::Matrix{Float64}
    r::Vector{Float64}
    z::Vector{Float64}
    phi::Vector{Float64}
    nr::Int64
    nz::Int64
    nphi::Int64
end

"""
rz_grid(rmin, rmax, nr, zmin, zmax, nz, phimin=0.0, phimax=0.0, nphi=1)

Creates interpolation grid.

#### Input arguments
## rmin - Minimum radius [cm]
## rmax - Maximum radius [cm]
## nr - Number of radii
## zmin - Minimum z value [cm]
## zmax - Maximum z value [cm]
## nz - Number of z values
## phimin - Minimum Phi value [rad]
## phimax - Maximum Phi value [rad]
## nphi - Number of Phi values 

#### Return Value
## Interpolation grid dictionary

####Example Usage
## ```julia
## julia> grid = rz_grid(0,200.0,200,-100,100,200,phimin=4*np.pi/3,phimax=5*np.pi/3,nphi=5)
## ```
"""
function rz_grid(rmin::Union{Int64,Float64}, rmax::Union{Int64,Float64}, nr::Int64, zmin::Union{Int64,Float64}, zmax::Union{Int64,Float64}, nz::Int64; phimin::Union{Int64,Float64}=0.0, phimax::Union{Int64,Float64}=0.0, nphi::Int64=1)
    dr = (rmax - rmin) / nr
    dz = (zmax - zmin) / nz
    dphi = (phimax - phimin) / nphi
    r = rmin .+ dr * (0:(nr-1))
    z = zmin .+ dz * (0:(nz-1))
    phi = phimin .+ dphi * (0:(nphi-1))

    r2d = repeat(r,1,nz)
    z2d = repeat(z,1,nr)'

    return grid(r2d,z2d,r,z,phi,nr,nz,nphi)
end

"""
read_ncdf(filepath; vars=nothing)

# Function description here
"""
function read_ncdf(filepath::String; wanted_keys=nothing)

    d = Dict()
    d["err"] = 1
    if isfile(filepath)
        d["err"] = 0
        NetCDF.open(filepath) do nc
            cdf_variables = nc.vars
            if !(wanted_keys==nothing)
                for wanted_key in wanted_keys
                    if wanted_key in keys(cdf_variables)
                        values = NetCDF.readvar(nc,wanted_key)
                        if () == size(values) # If the size of values is 0 (i.e. values is a scalar)
                            d[wanted_key] = first(values) # Take the 'first' element, parse a float and store it in d with key 'wanted_key'
                        else
                            if typeof(values) <: Vector{NetCDF.ASCIIChar}
                                values = reduce(*,map(x-> "$(x)",values)) # Concatenate all ASCII characters into one string
                            end
                            d[wanted_key] = values # Parse all elements in values as floats and store them as an array in d with key 'wanted_key'
                        end
                    end
                end
            else
                for (key,_) in cdf_variables
                    values = NetCDF.readvar(nc,key)
                    if () == size(values) # If the size of values is 0 (i.e values is a scalar)
                        d[key] = first(values) # Take the 'first' element, parse a float and store it in d with key 'wanted_key'
                    else
                        if typeof(values) <: Vector{NetCDF.ASCIIChar}
                            values = reduce(*,map(x-> "$(x)",values)) # Concatenate all ASCII characters into one string
                        end
                        d[key] = values # Parse all elements in values as floats and store them as an array in d with key 'wanted_key'
                    end
                end
            end
        end
    else
        error("FILE DOES NOT EXIST: "*filepath)
    end
    return d
end

"""
read_nubeam(filepath, mygrid)
- with example of non-default keyword arguments:
read_nubeam(filepath, mygrid; e_range=(5.0,500.0), p_range=(-0.75,0.75), btipsign=1, species=2, verbose=true, plotting=true)

Read the TRANSP NUBEAM-generated data in the .cdf-file 'filepath'. Use 'mygrid' to interpolate the fast-ion data onto.
NUBEAM-generated fast-ion distributions are given on an irrregular grid spiraling outwards from the magnetic axis. read_nubeam() will 
interpolate that data onto a regular (E,p,R,z) grid.

- If 'e_range' and/or 'p_range' are not specified, use TRANSP-data ranges by default. e_range should be specified in keV. If specified, only use TRANSP-data 
within the 'e_range'/'p_range' ranges to create the fast-ion distribution f(E,p,R,z).
- If 'btipsign' is not specified, assume -1 by default. 'btipsign' is the sign of the dot-product between the magnetic field and the plasma current.
- If 'species' is not specified, assume 1 by default. 'species' is the integer index of the fast-ion species in the TRANSP .cdf-file. Usually the first
fast-ion species '1' is wanted.
- If 'verbose' is not specified, assume false by default. If 'verbose' is set to true, read_nubeam() will talk a lot!
- If 'plotting' is not specified, assume false by default. If 'plotting' is set to true, read_nubeam() will plot spiral grid data and some test plots,
to enable verification that read_nubeam() is working correctly.
- If 'vverbose' is not specified, assume false by default. If 'vverbose' it set to true, read_nubeam() will talk a lot more (including printing for-loop information etc.)!
- If 'save_plots' is not specified, assume false by default. If 'save_plots' is set to true (and 'plotting' is set to true), the safety-check plots will be saved in .png format in your working directory.
Original function written in Python by Luke Stagner as part of FIDASIM.
"""
function read_nubeam(filepath::String, mygrid::grid; e_range::Union{Nothing,Tuple{Float64,Float64}}=nothing, p_range::Union{Nothing,Tuple{Float64,Float64}}=nothing, btipsign::Int64=-1, species::Int64=1, verbose::Bool = false, plotting::Bool = false, vverbose::Bool = false, save_plots::Bool = false)
    
    verbose && println("Loading data from TRANSP .cdf file... ")
    species_var = "SPECIES_$(species)" # The fast-ion species identification pattern for TRANSP data
    sstr = read_ncdf(filepath_distr, wanted_keys=[species_var])[species_var] # Will return a string with the fast-ion species of the TRANSP shot
    TRANSP_data = read_ncdf(filepath_distr, wanted_keys=["TIME","R2D","Z2D","E_"*sstr,"A_"*sstr,"F_"*sstr,"RSURF","ZSURF","BMVOL"]) # Will return a struct with the requested TRANSP data

    ngrid = length(TRANSP_data["R2D"]) # The number of (spiral) grid points on which the TRANSP data is given

    if length(TRANSP_data["TIME"])>1
        TRANSP_timepoint = TRANSP_data["TIME"][1] # The tokamak shot timepoint of the TRANSP-data. In seconds
    else
        TRANSP_timepoint = TRANSP_data["TIME"] # The tokamak shot timepoint of the TRANSP-data. In seconds
    end

    TRANSP_R_vector = TRANSP_data["R2D"] # 1D array. R values in cm 
    TRANSP_z_vector = TRANSP_data["Z2D"] # 1D array. z values in cm
    TRANSP_R_mesh = TRANSP_data["RSURF"] # 2D array. R values in cm
    TRANSP_z_mesh = TRANSP_data["ZSURF"] # 2D array. z values in cm
    bmvol = TRANSP_data["BMVOL"] # Volume of elements at (spiral) grid points. In cm^-3
    p_vector = TRANSP_data["A_"*sstr] # 1D array
    E_vector = TRANSP_data["E_"*sstr]*(1.0e-3) # 1D_array. Energy values in keV (after multiplication with (1.0e-3) factor)
    F_Epbm = (TRANSP_data["F_"*sstr])*(1.0e3) # 3D array. size(F_Epbm)=(length(E_vector), length(p_vector), length(bmvol)). Distribution values in (keV)^-1 cm^-3
    F_Epbm = map(x-> x > 0.0 ? 0.5*x : 0.0, F_Epbm) #Multiply elements by 0.5 to convert to pitch instead of solid angle d_omega/4pi

    if btipsign < 0
        verbose && println("Flipping fast-ion distribution in pitch since sign(J ⋅ B)< 0... ")
        F_Epbm = reverse(F_Epbm,dims=2) # Reverse all elements in dimension of pitch, to account for plasma current and magnetic field pointing in different directions
    end
    if e_range==nothing # If an energy range has not been specified as input with the keyword arguments...
        e_range = (minimum(E_vector), maximum(E_vector))
    end
    if p_range==nothing # If a pitch range has not been specified as input with the keyword arguments...
        p_range = (minimum(p_vector), maximum(p_vector))
    end

    we = findall(x-> x >= e_range[1] && x <= e_range[2], E_vector) # Find the indices of all energies within specified energy range e_range
    wp = findall(x-> x >= p_range[1] && x <= p_range[2], p_vector) # Find the indices of all pitches within specified pitch range p_range
    E_vector = E_vector[we] # Trim energy vector accordingly
    nenergy = length(E_vector)
    p_vector = p_vector[wp] # Trim pitch vector accordingly
    npitch = length(p_vector)
    F_Epbm = F_Epbm[we,:,:] # Trim fast-ion distribution accordingly
    F_Epbm = F_Epbm[:,wp,:]
    dE = abs(E_vector[2] - E_vector[1]) # Assume uniform energy grid...
    dp = abs(p_vector[2] - p_vector[1]) # Assume uniform pitch grid
    emin, emax = clamp(minimum(E_vector) - 0.5*dE, 0.0, minimum(E_vector)), maximum(E_vector) + 0.5*dE # Make sure that lower energy bound is larger than zero
    pmin, pmax = clamp(minimum(p_vector) - 0.5*dp, -1.0, minimum(p_vector)), clamp(maximum(p_vector)+0.5*dp, maximum(p_vector), 1.0) # Make sure pitch boundaries are within [-1.0, 1.0]

    verbose && println("Energy min/max: $((emin, emax))")
    verbose && println("Pitch min/max: $((pmin, pmax))")

    # Query (R,z) grid data on which to interpolate on
    query_nR = mygrid.nr # Number of query R points
    query_nz = mygrid.nz
    query_R = mygrid.r
    query_z = mygrid.z
    query_R_mesh = mygrid.r2d
    query_z_mesh = mygrid.z2d
    dR = abs(query_R[2] - query_R[1])
    dz = abs(query_z[2] - query_z[1])

    F_bm = dropdims(sum(F_Epbm,dims=(1,2))*dE*dp,dims=(1,2)) # The fast-ion distribution as a function of the bm volume elements (units: cm^-3)
    ntot = sum(F_bm .*bmvol) # Multiply F_bm with the bm volume elements to get the total number of fast ions. Units: (cm^-3)*(cm^3) = #fast ions
    verbose && println("Ntotal in phase space: $(ntot)")

    verbose && println("Creating the Delaunay tessellation... ")
    tri = DelaunayTessellation(length(TRANSP_R_vector)) # Instantiate a Delaunay tessellation of approximately length(TRANSP_R_vector)
    mnr = minimum([minimum(TRANSP_R_vector),minimum(query_R)]) # The definite minimum of every possible R coordinate in the algorithm
    mxr = maximum([maximum(TRANSP_R_vector),maximum(query_R)]) # The definite maximum of every possible R coordinate in the algorithm
    mnz = minimum([minimum(TRANSP_z_vector),minimum(query_z)]) # The definite minimum of every possible z coordinate in the algorithm
    mxz = maximum([maximum(TRANSP_z_vector),maximum(query_z)]) # The definite maximum of every possible z coordinate in the algorithm
    r_tri(x) = min_coord + (max_coord-min_coord)*(x-mnr)/(mxr-mnr) # DelaunayTesselation expects all points to be within [min_coord,max_coord] (given by VoronoiDelaunay.min_coord etc). We therefore need a conversion function.
    z_tri(y) = min_coord + (max_coord-min_coord)*(y-mnz)/(mxz-mnz) # DelaunayTesselation expects all points to be within [min_coord,max_coord] (given by VoronoiDelaunay.min_coord etc). We therefore need a conversion function.
    r_tri_inv(x_tri) = mnr + (mxr-mnr)*(x_tri-min_coord)/(max_coord-min_coord) # We need a function to convert back to normal R values
    z_tri_inv(y_tri) = mnz + (mxz-mnz)*(y_tri-min_coord)/(max_coord-min_coord) # We need a function to convert back to normal z values

    TRANSP_spiral_rz_iterator_tri = zip(r_tri.(TRANSP_R_vector), z_tri.(TRANSP_z_vector)) # Create an iterator with all the TRANSP (R,z) coordinates as tuples
    a = Point2D[Point(rr_tri, zz_tri) for (rr_tri,zz_tri) in TRANSP_spiral_rz_iterator_tri] # Put them all in a Point() object within a Point2D array
    push!(tri, a) # Feed all points to the Delaunay tessellation

    verbose && println("Mapping all triangles to original vertex indices... ")
    vertices = Dict{String,Tuple{Int64,Int64,Int64}}() # A hashmap. The hashes are triangle IDs and the values are tuples with vertices indices in the same order as the a,b and c fields of the DelaunayTriangle object returned by locate(tri, point). The DelaunayTesselation struct provided no way to keep track of what indices of the TRANSP_R_vector, TRANSP_z_vector inputs that make up the vertices of the triangles in the tessellation... So I had to do it myself!
    count = 1
    for t in tri # For triangles in the tessellation
        vverbose && println("Triangle $(count)")
        ind_a = findall(x-> x==(getx(geta(t)),gety(geta(t))), collect(TRANSP_spiral_rz_iterator_tri)) # Find the (R,z) data index of the first vertex to triangle t
        ind_b = findall(x-> x==(getx(getb(t)),gety(getb(t))), collect(TRANSP_spiral_rz_iterator_tri)) # Find the (R,z) data index of the second vertex to triangle t
        ind_c = findall(x-> x==(getx(getc(t)),gety(getc(t))), collect(TRANSP_spiral_rz_iterator_tri)) # Find the (R,z) data index of the third vertex to triangle t

        vverbose && println("ind_a: $(ind_a)    ind_b: $(ind_b)     ind_c: $(ind_c)")

        triangle_ID = "$(t._neighbour_a)_$(t._neighbour_b)_$(t._neighbour_c)" # No other triangle in the tessellation will have this ID. Thus making it easily look-up-able later
        vertices[triangle_ID] = (first(ind_a), first(ind_b), first(ind_c))
        vverbose && println("Triangle ID: "*triangle_ID)
        vverbose && println()
        count += 1
    end

    if plotting
        verbose && println("Plotting test plot of the Delaunay tessellation (requires Plots.jl package)... ")
        # Visualize the Delaunay tessellation (tessellation=the action of spanning an N-dimensional space with small building blocks. In this case: simplices. Simplice=an object in N-dimensional space with the smallest number of sides. In 2D, it's a triangel. In 3D, it's a tetrahedron. Etc.)
        x, y = getplotxy(delaunayedges(tri))
        test_plt_1 = Plots.plot(r_tri_inv.(x),z_tri_inv.(y); aspect_ratio=:equal, label="Delaunay tessellation")

        # Find a test point in the tessellation. Plot it, and the vertices of the surrounding triangle
        myp_R = mean(TRANSP_R_vector)
        myp_z = mean(TRANSP_z_vector)
        t = locate(tri, Point(r_tri(myp_R), z_tri(myp_z))) # Try to locate the query point inside of the tessellation
        verbose && println("-- Including test point $((myp_R,myp_z)) within the tessellation... ")
        test_plt_1 = Plots.scatter!([myp_R],[myp_z],label="Test point inside tessellation",ms=3.0)
        test_plt_1 = Plots.scatter!([r_tri_inv(getx(geta(t)))],[z_tri_inv(gety(geta(t)))],label="", ms=3.0)
        test_plt_1 = Plots.scatter!([r_tri_inv(getx(getb(t)))],[z_tri_inv(gety(getb(t)))],label="", ms=3.0)
        test_plt_1 = Plots.scatter!([r_tri_inv(getx(getc(t)))],[z_tri_inv(gety(getc(t)))],label="", ms=3.0)

        # Test hashmap on test point
        tID = "$(t._neighbour_a)_$(t._neighbour_b)_$(t._neighbour_c)"
        ia, ib, ic = vertices[tID]
        x = [r_tri(myp_R), z_tri(myp_z)] # Query point in tri coordinates
        xa = [getx(geta(t)), gety(geta(t))] # Vertex a of triangle t
        xb = [getx(getb(t)), gety(getb(t))] # Vertex b of triangle t
        xc = [getx(getc(t)), gety(getc(t))] # Vertex c of triangle t

        μ = [(1/sum(abs2,x-xa)),(1/sum(abs2,x-xb)),(1/sum(abs2,x-xc))] 
        μ = (1/sum(μ)) .* μ # Barycentric weights. sum(μ)=1 must hold
        fdens_test = μ[1]*F_bm[ia]+μ[2]*F_bm[ib]+μ[3]*F_bm[ic]
        verbose && println("-- The barycentric weights, and corresponding fast-ion distribution values, given the test point $((myp_R,myp_z)) are...")
        verbose && println("---- μ[1], F_bm[a]: $((μ[1],F_bm[ia]))")
        verbose && println("---- μ[2], F_bm[b]: $((μ[2],F_bm[ib]))")
        verbose && println("---- μ[3], F_bm[c]: $((μ[3],F_bm[ic]))")
        verbose && println("---- Interpolated value of F_bm at $((myp_R,myp_z)) given the barycentric weights: $(fdens_test)") # The interpolation makes sense! I think.
    end

    # Test hashmap on test point
    verbose && println("Creating nearest neighbour (NN) tree for extrapolation of points outside the tessellation... ")
    brutetree = BruteTree(hcat(TRANSP_R_vector,TRANSP_z_vector)') # There are also other tree types. However, BruteTree was deemed good enough
    if plotting
        # Find a test point outside of the tessellation. Plot it, including its NN
        myNN_R = maximum(TRANSP_R_vector) + 0.05*abs(maximum(TRANSP_R_vector)-minimum(TRANSP_R_vector))
        myNN_z = mean(TRANSP_z_vector)+0.05*abs(maximum(TRANSP_z_vector)-minimum(TRANSP_z_vector))
        verbose && println("-- Including test point $((myNN_R, myNN_z)) outside of the tessellation... ")
        idx, dist = nn(brutetree, [myNN_R, myNN_z]) # nn=nearest neighbour. brutetree is just a tree to brute-force search for the nn to the test point (myNN_R, myNN_z)
        test_plt_1 = Plots.scatter!(TRANSP_R_vector, TRANSP_z_vector; aspect_ratio=:equal, label="")
        test_plt_1 = Plots.scatter!([myNN_R],[myNN_z], label="Test point outside of tessellation (qp)")
        test_plt_1 = Plots.scatter!([TRANSP_R_vector[idx]], [TRANSP_z_vector[idx]], label="NN to qp", xlabel="R [cm]", ylabel="z [cm]")
        display(test_plt_1) # Display test plot 1
        save_plots && png("test_plt_1")
        save_plots && println("Plot saved!")
    end

    F_Rz = zeros(query_nR,query_nz)
    F_EpRz = zeros(nenergy,npitch,query_nR,query_nz)
    count = 1
    verbose && println("Interpolating fast-ion distribution onto all query points... ")
    for (ir,rr) in enumerate(query_R), (iz,zz) in enumerate(query_z)
        vverbose && println("$(count)/$(length(query_R)*length(query_z)): (R,z)=($(rr),$(zz))")
        t = locate(tri, Point(r_tri(rr), z_tri(zz))) # Try to locate the query point inside of the tessellation

        if isexternal(t) == true # If the query point was outside of the tessellation...
            idx, dist = nn(brutetree, [rr, zz]) # Use brute-force nearest-neighbour extrapolation instead! Notice that here we don't need to r_tri() and z_tri() functions
            F_Rz[ir,iz] = F_bm[idx]
            F_EpRz[:,:,ir,iz] = F_Epbm[:,:,idx]
        else
            tID = "$(t._neighbour_a)_$(t._neighbour_b)_$(t._neighbour_c)" # Contruct the triangle ID from the order of the neighbours
            ia, ib, ic = vertices[tID] # Look up the original indices of the vertices
            x = [r_tri(rr), z_tri(zz)] # Query point in tri coordinates
            xa = [getx(geta(t)), gety(geta(t))] # Vertex a of triangle t
            xb = [getx(getb(t)), gety(getb(t))] # Vertex b of triangle t
            xc = [getx(getc(t)), gety(getc(t))] # Vertex c of triangle t

            μ = [(1/sum(abs2,x-xa)),(1/sum(abs2,x-xb)),(1/sum(abs2,x-xc))] 
            μ = (1/sum(μ)) .* μ # Barycentric weights. sum(μ)=1 must hold

            # Barycentric interpolation. All components of μ sum to one. 'Barycentric' simply means 'center-of-mass-like'. https://en.wikipedia.org/wiki/Barycentric_coordinate_system 
            F_Rz[ir,iz] = μ[1]*F_bm[ia]+μ[2]*F_bm[ib]+μ[3]*F_bm[ic]
            F_EpRz[:,:,ir,iz] = μ[1]*F_Epbm[:,:,ia]+μ[2]*F_Epbm[:,:,ib]+μ[3]*F_Epbm[:,:,ic]
        end
        count += 1
    end 

    F_Rz = map(x-> x<0.0 ? 0.0 : x, F_Rz) # Any negative values? Set them to zero instead!

    # Correct for points outside of seperatrix. They should all be zero
    verbose && println("Identifying query points outside of the separatrix... ")
    R_maxis = mean(TRANSP_R_mesh[:,1]) # The R coordinate of the magnetic axis
    z_maxis = mean(TRANSP_z_mesh[:,1])
    R_sep = TRANSP_R_mesh[:,end] # The R coordintes of the separatrix
    z_sep = TRANSP_z_mesh[:,end]

    
    if plotting # Double-check axis and sep values via plotting
        verbose && println("-- Plotting a safety check plot for data points, mag. axis and separatrix... ")
        x, y = getplotxy(delaunayedges(tri))
        test_plt_2 = Plots.plot(r_tri_inv.(x),z_tri_inv.(y); aspect_ratio=:equal, label="Delaunay tessellation")
        test_plt_2 = Plots.scatter!(TRANSP_R_vector,TRANSP_z_vector,label="Nodes")
        test_plt_2 = Plots.scatter!([R_maxis],[z_maxis],label="Mag. axis")
        test_plt_2 = Plots.scatter!(R_sep,z_sep,label="Separatrix", xlabel="R [cm]", ylabel="z [cm]")
        display(test_plt_2) # Display test plot 2
        save_plots && png("test_plt_2")
        save_plots && println("Plot saved!")
    end

    R_bdry = R_sep .- R_maxis # The R coordinates of the separatrix, relative to the magnetic axis
    z_bdry = z_sep .- z_maxis

    r_bdry = sqrt.(R_bdry.^2 .+ z_bdry.^2) # The r coordinates (local toroidal coordinate system. http://fusionwiki.ciemat.es/wiki/Toroidal_coordinates)
    theta_bdry = atan.(z_bdry,R_bdry) # The θ coordinates from arctan of (z-z_axis)/(R-R_axis) values
    theta_bdry = map(x-> x < 0.0 ? (x + 2*pi) : x,theta_bdry) # Ensure θ coordinates within [0,2π]
    w = sortperm(theta_bdry) # Sort theta_bdry but return only the indices that would give the sorted order
    theta_bdry = theta_bdry[w]
    r_bdry = r_bdry[w]
    w = unique(i -> theta_bdry[i], 1:length(theta_bdry)) # Keep only unique elements
    theta_bdry = theta_bdry[w]
    r_bdry = r_bdry[w]

    itp = Interpolations.interpolate((theta_bdry,), r_bdry, Gridded(Linear())) # Create an interpolation object, to be able to interpolate new values of r from input theta values
    etpf = extrapolate(itp, Flat()) # If outside of domain, just return flat (constant) values (flat=last known data value before extrapolation domain)


    R_query_pts = query_R_mesh .- R_maxis # Same procedure as for data points, but for query points
    z_query_pts = query_z_mesh .- z_maxis
    r_query_pts = sqrt.(R_query_pts.^2 .+ z_query_pts.^2)
    theta_query_pts = atan.(z_query_pts,R_query_pts)
    theta_query_pts = map(x-> x < 0.0 ? (x + 2*pi) : x, theta_query_pts) # [0,2pi]

    # For all query theta angles, find the cooresponding r points on the separatrix
    r_query_bdry = zeros(size(theta_query_pts))
    for (i,theta_query_pt) in enumerate(theta_query_pts)
        r_query_bdry[i] = etpf(theta_query_pt)
    end

    w = findall(x-> x >= 0.0, (r_query_pts - r_query_bdry) .- 2) # If the query points are more than 2 cm outside of the separatrix...
    F_Rz[w] .= 0.0 # ...make sure that the fast-ion distribution is zero there and outwards (> 2cm outside of the separatrix)
    F_EpRz[:,:,w] .= 0.0 # ...make sure that the fast-ion distribution is zero there and outwards (> 2cm outside of the separatrix)

    # Enforce correct normalization (sum(distribution) = #total number of fast ions) for interpolated fast-ion distributions
    verbose && println("Enforcing correct normalization for interpolated fast-ion distribution (sum(f)=#total number of fast ions)... ")
    ntot_denf = 2*pi*dR*dz*sum(reshape(query_R,(length(query_R),1)) .*F_Rz)
    F_Rz = (ntot/ntot_denf) .* F_Rz 
    ntot_fbm = (2*pi*dE*dp*dR*dz)*sum(reshape(query_R,(1,1,length(query_R),1)) .*F_EpRz)
    F_EpRz = (ntot/ntot_fbm) .*F_EpRz

    verbose && println("Collecting all data in dictionary and returning it as output... ")
    fbm_dict = Dict()
    fbm_dict["species"] = sstr # Fast-ion species data
    fbm_dict["time"] = TRANSP_timepoint # The tokamak shot timepoint in seconds
    fbm_dict["ntot"] = ntot # Total number of fast ions in distribution
    fbm_dict["energy"] = E_vector # Fast-ion energies in keV
    fbm_dict["pitch"] = p_vector # Fast-ion pitches
    fbm_dict["R"] = query_R # Fast-ion R values in centimeters
    fbm_dict["z"] = query_z # Fast-ion z values in centimeters
    fbm_dict["F_EpRz"] = F_EpRz # Fast-ion distribution as 4D array, in (keV)^-1 (cm^-3)
    fbm_dict["F_Rz"] = F_Rz # Fast-ion distribution as 2D array in (R,z), in cm^-3
    fbm_dict["R_sep"] = R_sep # R coordinates of separatrix
    fbm_dict["z_sep"] = z_sep # z coordinates of separatrix
    fbm_dict["btipsign"] = btipsign
    fbm_dict["data_source"] = filepath

    return fbm_dict
end

############################################################# Execute the script
# Change directory to the OWCF extra/ folder

# Split filepath_distr by using '/' as delimeter
verbose && println("Examining filepath_distr... ")
fd_parts = split(filepath_distr,'/')

# Assume last element is the actual filename
filename_distr = fd_parts[end]

# Remove filename extension
filename_distr = (split(filename_distr,'.'))[1]

# Define the (R,z) grid
verbose && println("Computing (R,z) grid for fast-ion distribution... ")
mygrid = rz_grid(Rmin, Rmax, nR, zmin, zmax, nz)

# Read the TRANSP computed fast-ion distribution, and convert it to (E,p,R,z) format
verbose && println("Converting fast-ion distribution from TRANSP .cdf-file spiral format to (E,p,R,z) format... ")
F_EpRz_dict = read_nubeam(filepath_distr,mygrid; e_range=e_range, p_range=p_range, species=species, btipsign=btipsign, verbose=verbose, plotting=plotting, vverbose=vverbose, save_plots=save_plots)

# Save the results in .h5 file
global filepath_output_orig = folderpath_out*filename_distr
global filepath_output = deepcopy(filepath_output_orig)
global count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output = filepath_output_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_output = filepath_output*".jld2"
verbose && println("Saving the results at "*filepath_output*"... ")
jldopen(filepath_output, true, true, false) do file
    write(file, "species", F_EpRz_dict["species"])
    write(file, "timepoint", F_EpRz_dict["time"])
    write(file, "ntot", F_EpRz_dict["ntot"])
    write(file, "F_EpRz", F_EpRz_dict["F_EpRz"])
    write(file, "F_Rz", F_EpRz_dict["F_Rz"])
    write(file, "energy", F_EpRz_dict["energy"])
    write(file, "pitch", F_EpRz_dict["pitch"])
    write(file, "R", F_EpRz_dict["R"])
    write(file, "z", F_EpRz_dict["z"])
    write(file, "R_sep", F_EpRz_dict["R_sep"])
    write(file, "z_sep", F_EpRz_dict["z_sep"])
    write(file, "signOfJdotB", F_EpRz_dict["btipsign"])
    write(file, "orig_TRANSP_filepath", F_EpRz_dict["data_source"])
end

verbose && println("~~~~~~ getEpRzFIdistrFromTRANSP.jl completed successfully! ~~~~~~")