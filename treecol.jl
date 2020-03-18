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
const NSIDE = 4
using Healpix

#include("OctTree.jl")
using OctTree
#include("ChemistryNetwork.jl")
using ChemistryNetwork



T = Float64
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

    return N_gas, pos_gas, vel_gas, rho, u, m_gas, hsml, id_gas, boxsize, time
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

const ANGLE = 0.3

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
    file_path = "./isocloud_N1e4/"
    #file_path = "./nH10_box150pc_S4_N1e6_myIC/"
    #file_path = "./nH1_box250pc_S2_N1e5_myIC/"
    snap = "013"
    filename = file_path * "/snap_" * snap * ".hdf5"
    Npart, pos, vel, rho, u, mass, hsml, id_gas, boxsize, time = read_snap(filename);
    X = vec2svec(pos);
    boxsizes = SVector(1.0, 1.0, 1.0) .* (boxsize)

    #X .*= BOXSIZE
    #push!(X,SVector(0.,0.,0.) .+ 0.01)


    #nH = 1000.
    #temp = 50.
    ξ = 1.3e-16 #H2
    IUV = 1.0
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

    abund_all = [zeros(ChemistryNetwork.N_spec) for _ in 1:Npart]

    for i in 1:Npart
        #abund = abund_all[:,1]
        init_abund(abund_all[i], Zp)
    end

    #hsml = ones(Npart) .* (0.5 * boxsizes[1])
    NH_part = zeros(Npart)
    NH2_part = zeros(Npart)
    NCO_part = zeros(Npart)
    #ga_out = TreeGather{T}()
    ga_out = Vector{TreeGather{Float64}}(undef,Npart)
    tree_out = nothing

    #ITER = 1
    Nstep = 10
    time_max = 1e8 #unit = yr
    dt = time_max / Nstep

    for j in 1:Nstep
        @. mass_H2 = mass * getindex(abund_all,iH2)
        @. mass_CO = mass * getindex(abund_all,iCO)
        #@show mass_H2[1], mass_CO[1]
        println("j = ", j)
        println("buildtree...")
        @time tree = buildtree(X, hsml, mass, mass_H2, mass_CO, boxsizes);
        println("done!")
        #@time for i in 1:Npart
        println("loop over particles...")
        #@time for i in 1:Npart
        @time @threads for i in 1:Npart
            #i%100 == 0 ? println("i = ", i) : nothing
            ga = TreeGather{T}()
            treewalk(ga,X[i],tree,ANGLE,boxsizes)
            #column_all = (ga.column_all);
            NH = median(ga.column_all) * fac_col
            NH2 = median(ga.column_H2) * fac_col
            NCO = median(ga.column_CO) * fac_col
            NC = 0.0
            #@show NH, NH2, NCO, NC
            #if i==Npart
            #    ga_out = ga
            #end            
            ga_out[i] = ga
            NH_part[i] = NH
            NH2_part[i] = NH2
            #NCO_part[i] = NCO
            NCO_part[i] = i
            solve_equilibrium_abundances(abund_all[i], dt, nH[i], temp[i], NH, NH2, NCO, NC, ξ, IUV, Zp)
        end
        println("done!")
        tree_out = tree

        # write to file
        open(file_path * "/chem3D-N" * string(Npart) * "-" * string(j) *".out","w") do f
            serialize(f, NH_part)
            serialize(f, NH2_part)
            serialize(f, NCO_part)
            serialize(f, abund_all)
            serialize(f, ga_out[Npart].column_all)
            serialize(f, X)
            serialize(f, dict)
        end
    end
    return abund_all, ga_out, X, NH_part, NH2_part, NCO_part, tree_out
end

abund_all, ga, X, NH, NH2, NCO, tree = solve_chem_all_particles();
0
#=
xmin=-2
xmax=3
ymin=1
ymax=7
cc,xb,yb,im=hist2D(log10.(rho.*404), log10.(u.*100), bins=[100,100], range=[[xmin, xmax],[ymin, ymax]])
clf()
imshow(transpose(log10.(cc)), origin="lower", extent=(xmin, xmax, ymin, ymax))
=#

#=
clf()
r = [sqrt(sum((X[i].- 0.5*0.01).^2)) for i in 1:Npart]
plot(r, NH, "o", ms=0.5, mfc="none", label="NH_tot", c="red")
plot(r, @.( fac_col*median(getfield(ga,:column_all)) ), "o", ms=0.5, mfc="none", label="NH_tot", c="blue")
gcf()
=#

#m = Map{Float64, RingOrder}(NSIDE);
#m.pixels[:] = log10.(ga.column_all .+ 1e-2*minimum(ga.column_all[ga.column_all.!=0.0]));
#m.pixels[:] = log10.(ga.column_all .* fac_col) ;
#Plots.plot(m, clim=(20.5,24))
#Plots.plot(m)

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


#=
fig = figure("", figsize=(12,8))
subplot(121)
plot(r,getindex.(abund_all,dict["C+"]),".",ms=1)
plot(r,getindex.(abund_all,dict["C"]) ,".",ms=1)
plot(r,getindex.(abund_all,dict["CO"]),".",ms=1)
xscale("log")
yscale("log")
xlabel("radius [pc]")
ylabel("x_i")

subplot(122)
plot(NH,getindex.(abund_all,dict["C+"]),".",ms=1)
plot(NH,getindex.(abund_all,dict["C"]) ,".",ms=1)
plot(NH,getindex.(abund_all,dict["CO"]),".",ms=1)
xscale("log")
yscale("log")
xlabel("NH [cm^-2]")
ylabel("x_i")

savefig("PDR.png")

=#

# write to file
#=
using Serialization
open("PDR3D.out","w") do f
    serialize(f, NH)
    serialize(f, NH2)
    serialize(f, NCO)
    serialize(f, abund_all)
    serialize(f, ga.column_all)
    serialize(f, X)
    serialize(f, dict)
end
=#

#===
p = boxsizes .* 0.5

#ga = TreeGather{T}(0,zeros(NSIDE^2*12),zeros(NSIDE^2*12),zeros(NSIDE^2*12),[],[])


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
