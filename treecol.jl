push!(LOAD_PATH, pwd())
using HDF5
using StaticArrays
using Statistics
#using PyPlot
using LinearAlgebra
using .Threads
using Serialization
using Random
#import StatsBase

#import Plots

#include("OctTree.jl")
using OctTree
include("ChemistryNetwork.jl")
#using ChemistryNetwork

using Healpix


const T = Float64
const N = 3
const GEO_FAC = T(4.0 / 3.0 * pi)
const KERNELCONST = T(16.0/pi)
@inline ramp(x) = max(0, x);
@inline function kernel_cubic(x::T) where {T}
    return T(KERNELCONST * (ramp(1.0 - x)^3 - 4.0 * ramp(0.5 - x)^3))
end

function read_snap(filename)

    T=Float64
    header::Dict = h5readattr(filename, "/Header")
    boxsize::T = header["BoxSize"]
    time::T    = header["Time"]

    N_gas::Int64 = header["NumPart_ThisFile"][1]

    pos_gas::Matrix{T} = h5read(filename, "PartType0/Coordinates");
    vel_gas::Matrix{T} = h5read(filename, "PartType0/Velocities");
    rho::Vector{T}     = h5read(filename, "PartType0/Density");
    u::Vector{T}       = h5read(filename, "PartType0/InternalEnergy");
    m_gas::Vector{T}   = h5read(filename, "PartType0/Masses");
    hsml::Vector{T}    = h5read(filename, "PartType0/SmoothingLength");
    #scal::Vector{T}       = h5read(filename, "PartType0/PassiveScalarField");

    id_gas::Vector{Int64} = h5read(filename, "PartType0/ParticleIDs");

    abund::Matrix{T} = h5read(filename, "PartType0/ChemicalAbundancesSG");
    fH2::Vector{T}    = h5read(filename, "PartType0/ShieldingFactorH2");
    fdust::Vector{T}    = h5read(filename, "PartType0/ShieldingFactorDust");

    col_tot::Matrix{T}    = h5read(filename, "PartType0/TreecolColumnDensitiesAll");
    col_H2::Matrix{T}    = h5read(filename, "PartType0/TreecolColumnDensitiesH2");
    col_CO::Matrix{T}    = h5read(filename, "PartType0/TreecolColumnDensitiesCO");

    Tdust::Vector{T}    = h5read(filename, "PartType0/DustTemperature");

    N_star::Int64 = header["NumPart_ThisFile"][5]

    pos_star::Matrix{T} = h5read(filename, "PartType4/Coordinates");
    vel_star::Matrix{T} = h5read(filename, "PartType4/Velocities");
    m_star::Vector{T}   = h5read(filename, "PartType4/Masses");
    sftime::Vector{T}      = h5read(filename, "PartType4/StellarFormationTime");
    id_star::Vector{Int64} = h5read(filename, "PartType4/ParticleIDs");

    return N_gas, pos_gas, vel_gas, rho, u, m_gas, hsml, id_gas,
        abund, fH2, fdust, col_tot, col_H2, col_CO, Tdust,
        N_star, pos_star, vel_star, m_star, sftime, id_star,
        boxsize, time
end

function vec2svec(vec::Matrix{T}) where {T}
    svec = [SVector{3,T}(vec[:,i]) for i in 1:size(vec,2)]
end
function mat2smat(mat::Array{T,3}) where {T}
    smat = [SMatrix{3,3,T}(mat[:,:,i]) for i in 1:size(mat,3)]
end




#const BOXSIZE_X = BOXSIZE_Y = BOXSIZE_Z = header["BoxSize"]
#const boxsizes = SVector{N,T}(BOXSIZE_X, BOXSIZE_Y, BOXSIZE_Z)



const N=3
#const Npart = 1000;

#const BOXSIZE = 1e-2
#const BOXSIZE = header["BoxSize"]

#double the boxsize for non-periodic BC
#this is a hack; should do it properly in OctTree in the future
#const boxsizes = SVector(1.0, 1.0, 1.0) .* (BOXSIZE*2.0)
#const boxsizes = SVector(1.0, 1.0, 1.0) .* (BOXSIZE)

const ANGLE = 0.7
const ShieldingLength = 0.1

const XH = 0.71
const BOLTZMANN=1.3806e-16
const PROTONMASS=1.6726e-24
const GRAVCON=6.67e-8
const UnitMass_in_g = 1.989e43
const UnitLength_in_cm    = 3.085678e21
const UnitTime_in_s = 3.08568e+16

const UnitDensity_in_cgs = UnitMass_in_g / UnitLength_in_cm^3
const UnitDensity_in_pccm = UnitDensity_in_cgs/PROTONMASS
const Year_in_s = 31556926.

const fac_col = (UnitMass_in_g/UnitLength_in_cm^2)*(XH/PROTONMASS)



function solve_chem_all_particles()
    #X = [SVector{N}( rand(MersenneTwister(i),N) ) for i in 1:Npart-1]
    #push!(X,SVector(0.,0.,0.))
    #push!(X,SVector(0.5,0.5,0.5))

    #file_path = "./"
    #file_path = "../../../simulations/isocloud_N1e4/"
    #snap = "013"
    file_path = "../../../simulations/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft4_SFLJ4_eff0p1_stoIMFfix"
    snap = "1000"
    fname = file_path * "/snap_" * snap * ".hdf5"
    Npart, pos, vel, rho, u, mass, hsml, id_gas,
        abund, fH2, fdust, col_tot, col_H2, col_CO, Tdust,
        N_star, pos_star, vel_star, m_star, sftime, id_star,
        boxsize, time = read_snap(fname);
    #println("boxsize = ", boxsize)

    boxsizes = SVector(1.0, 1.0, 20.0)
    center = SVector(0.5, 0.5, 0.0)
    topnode_length = SVector(10., 10., 10.)


    #down-sampling
    #=
    Ns::Int = 100
    rng = MersenneTwister(1234);
    ridx = randperm(rng, Npart);
    Npart = div(Npart,Ns)
    ridx = ridx[1:Npart]
    pos=pos[:,ridx]
    vel=vel[:,ridx]
    rho=rho[ridx]
    u=u[ridx]
    mass=mass[ridx]
    hsml=hsml[ridx]
    id_gas=id_gas[ridx]
    abund=abund[:,ridx]
    mass .*= Ns
    =#

    X = vec2svec(pos);
    #X .*= BOXSIZE
    #push!(X,SVector(0.,0.,0.) .+ 0.01)

    dt = 3e-2   #30Myr
    idx = @. ( time - sftime < dt )
    SFR = 1e10 * sum(abs.(m_star[idx])) / (dt*1e9)
    facSFR = SFR / 5e-3
    @show facSFR, SFR

    ξ = 1.3e-16 * facSFR #H2
    IUV = 1.0 * facSFR
    Zp = 1.0

    #=
    rho = nH * PROTONMASS / XH  / UnitDensity_in_cgs #code units

    Vbox = BOXSIZE^3
    Mbox = Vbox * rho #code unit

    mass = ones(Npart) .* (Mbox / Npart)
    =#

    nH = rho .* (XH * UnitDensity_in_cgs / PROTONMASS)
    temp = u ./ (1.5 * BOLTZMANN /PROTONMASS / 1e10)

    mass_H2 = zeros(Npart)
    mass_CO = zeros(Npart)

    #NH = NH2 = NCO = NC = 0.0 #don't declear here!! otherwise they won't be thread-safe

    abund_all = [zeros(N_spec) for _ in 1:Npart]


    for i in 1:Npart
        #abund = abund_all[:,1]
        #xneq = SVector{1,T}([abund[1,i]])
        xneq = SVector{N_neq,T}([abund[1,i], abund[2,i]])
        #use C+ & CO from simulations as the initial guess
        #abund_all[i][dict["CO"]] = abund[3,i]
        #abund_all[i][dict["C+"]] = abC_s * Zp - abund_all[i][dict["CO"]]
        #make sure the abC is consistent between simulations & postprocessing
        abund_all[i][dict["C+"]] = abC_s * Zp
        #initial guess for H2 & H+ (useful for calcultating steady state H2)
        abund_all[i][dict["H2"]] = abund[1,i]
        abund_all[i][dict["H+"]] = abund[2,i]
        init_abund(abund_all[i], Zp, xneq)
    end

    #hsml = ones(Npart) .* (0.5 * boxsizes[1])
    NH_eff = zeros(Npart)
    NH2_eff = zeros(Npart)
    NCO_eff = zeros(Npart)
    #ga_out = TreeGather{T}()
    ga_out = Vector{TreeGather{Float64}}(undef,Npart)
    tree_out = nothing

    #ITER = 1
    Nstep = 4
    time_max = 1e9 #unit = yr
    dt = time_max / Nstep
    NCpix = zeros(NPIX)

    for j in 1:Nstep
        @. mass_H2 = mass * getindex(abund_all,iH2)
        @. mass_CO = mass * getindex(abund_all,iCO)
        #@show mass_H2[1], mass_CO[1]
        println("j = ", j)
        println("buildtree...")
        #@time tree = buildtree(X, hsml, mass, mass_H2, mass_CO, boxsizes);
        @time tree = buildtree(X, hsml, mass, mass_H2, mass_CO, center, topnode_length);
        println("done!")
        #@time for i in 1:Npart
        println("loop over particles...")
        #@time for i in 1:Npart
        @time @threads for i in 1:Npart
            if nH[i] < 1.0 || temp[i] > 3e3 continue end
            if i%100 == 0 print("  i=", i) end
            #i%100 == 0 ? println("i = ", i) : nothing
            ga = TreeGather{T}()
            treewalk(ga,X[i],tree,ANGLE,ShieldingLength,boxsizes)
            #column_all = (ga.column_all);
            ga_out[i] = ga
            NH = ga.column_all .* fac_col
            NH2 = ga.column_H2 .* fac_col
            NCO = ga.column_CO .* fac_col
            NC = 0.0
            facNHtoAV = 5.35e-22
            #NH_eff[i] = median(NH)
            NH_eff[i] = -log(mean(exp.(-NH.*facNHtoAV))) / facNHtoAV
            NH2_eff[i] = -log(mean(exp.(-NH2.*facNHtoAV))) / facNHtoAV
            NCO_eff[i] = -log(mean(exp.(-NCO.*facNHtoAV))) / facNHtoAV
            #if median(NH) > 0
            #    @show ga.column_all, ga.column_H2, ga.column_CO
            #end
            #NH = ga.column_all .* fac_col
            #NH2 = ga.column_H2 .* fac_col
            #NCO = ga.column_CO .* fac_col
            xneq = SVector{N_neq,T}([abund[1,i], abund[2,i]])

            par = Par{NPIX,T}(nH[i], temp[i], ξ, IUV, Zp,
                SVector{NPIX,T}(NH),
                SVector{NPIX,T}(NH2),
                SVector{NPIX,T}(NCO),
                SVector{NPIX,T}(NCpix), xneq)
            #dt = 1e9 / (Zp * nH[i]) #in years
            solve_equilibrium_abundances(abund_all[i], dt, par)
        end
        println("done!")
        tree_out = tree


        # write to file
        abund_all_arr = zeros(N_spec,Npart)
        for i in eachindex(abund_all)
            abund_all_arr[:,i] = abund_all[i]
        end
        #T = Float64
        fname = file_path * "/chem3D-neqH2Hp-" * string(j) *".hdf5"
        fid=h5open(fname,"w")
        grp_head = g_create(fid,"Header");
        attrs(fid["Header"])["all_species"] = all_species
        attrs(fid["Header"])["time"] = time
        attrs(fid["Header"])["facSFR"] = facSFR
        grp_part = g_create(fid,"Chemistry");
        h5write(fname, "Chemistry/Abundances"     , abund_all_arr)
        h5write(fname, "Chemistry/ID"             , id_gas)
        h5write(fname, "Chemistry/NH_eff"         , NH_eff)
        h5write(fname, "Chemistry/NH2_eff"        , NH2_eff)
        h5write(fname, "Chemistry/NCO_eff"        , NCO_eff)
        close(fid)

        open(file_path * "/chem3D-neqH2Hp-" * string(j) *".out","w") do f
            serialize(f, NH_eff)
            serialize(f, NH2_eff)
            serialize(f, NCO_eff)
            serialize(f, abund_all)
            #serialize(f, ga_out[Npart].column_all)
            serialize(f, rho)
            serialize(f, u)
            serialize(f, dict)
        end

    end
    return abund_all, ga_out, X, NH_eff, NH2_eff, NCO_eff, tree_out
end

abund_all, ga, X, NH, NH2, NCO, tree = solve_chem_all_particles();
0



#=
ix,iy = 1,2
clf()
fig = figure("",figsize=(8,8))
#plot_quadtree(tree, ix, iy)
scatter(getindex.(X,ix), getindex.(X,iy), c="blue", marker="o", s=0.3)
plot_treewalk(ga,ix,iy)
p = X[Npart]
plot(p[ix],p[iy],"o")
gcf()
#savefig("particles.png")
=#


#===
p = boxsizes .* 0.5

m = Map{Float64, RingOrder}(NSIDE);
m.pixels[:] = log10.(column_all)
#Plots.plot(m, clim=(3e4,8e4))
Plots.plot(m)

column_true = zeros(12*NSIDE^2)
@time for i in 1:length(X)
    dx = OctTree.nearest.(X[i] - p, boxsizes)
    ipix = OctTree.vec2pix(dx)
    #ga.column_all[ipix] += node.mass
    area = 4*pi/(NSIDE^2*12.) * sum(dx.^2)
    column_true[ipix] += mass[i] / area
end

m.pixels[:] = log10.(column_true)
Plots.plot(m)
===#
