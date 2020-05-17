using Serialization
using StaticArrays
using PyPlot
using HDF5
using Random

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
    fH2ss::Vector{T}    = h5read(filename, "PartType0/ShieldingFactorH2");
    fH2dust::Vector{T}    = h5read(filename, "PartType0/ShieldingFactorDust");

    col_tot::Matrix{T}    = h5read(filename, "PartType0/TreecolColumnDensitiesAll");
    col_H2::Matrix{T}    = h5read(filename, "PartType0/TreecolColumnDensitiesH2");
    col_CO::Matrix{T}    = h5read(filename, "PartType0/TreecolColumnDensitiesCO");

    Tdust::Vector{T}    = h5read(filename, "PartType0/DustTemperature");

    return N_gas, pos_gas, vel_gas, rho, u, m_gas, hsml, id_gas,
        abund, fH2ss, fH2dust, col_tot, col_H2, col_CO, Tdust, boxsize, time
end

num = "1000"
file_path = "../../../simulations/isocloud_N1e4/"
#file_path = "../../../simulations/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft4_SFLJ4_eff0p1_stoIMFfix"

#Ngas, pos, vel, rho, u, mass, hsml, id,
#abund, fH2ss, fH2dust, col_tot, col_H2, col_CO, Tdust, boxsize, time =
#        read_snap(file_path * "/snap_" * num * ".hdf5");

# load from file
fname = file_path * "/chem3D-neqH2Hp-1"

open(fname*".out","r") do f

    global NH, NH2, NCO, abund_all, column_all, rho, u, d
    NH = deserialize(f)
    NH2 = deserialize(f)
    NCO = deserialize(f)
    abund_all = deserialize(f)
    #column_all = deserialize(f)
    rho = deserialize(f)
    u = deserialize(f)
    d = deserialize(f)
end

nH = rho .* 287

#boxsize = 0.15
#r = [sqrt(sum((X[i].- 0.5*boxsize).^2)) for i in eachindex(X)]

ms = 0.3
clf()

fig, ax = PyPlot.subplots(2, 2, figsize=(10,8))

ax[1,1].plot(nH,getindex.(abund_all,d["H"]) , "o", mfc="none", ms=ms, label="H")
ax[1,1].plot(nH,getindex.(abund_all,d["H2"]), "o", mfc="none", ms=ms, label="H2")
ax[1,1].set_xscale("log")
ax[1,1].set_yscale("log")
ax[1,1].set_xlabel("nH [cm^-3]")
ax[1,1].set_ylabel("x_i")
ax[1,1].set_ylim(1e-3,2e0)
ax[1,1].legend(loc="best", frameon=false, markerscale=5.0/ms)

ax[1,2].plot(NH,getindex.(abund_all,d["H"]) , "o", mfc="none", ms=ms, label="H")
ax[1,2].plot(NH,getindex.(abund_all,d["H2"]), "o", mfc="none", ms=ms, label="H2")
ax[1,2].set_xscale("log")
ax[1,2].set_yscale("log")
ax[1,2].set_xlabel("NH_tot  [cm^-2]")
ax[1,2].set_ylabel("x_i")
ax[1,2].set_ylim(1e-3,2e0)
ax[1,2].legend(loc="best", frameon=false, markerscale=5.0/ms)

ax[2,1].plot(nH,getindex.(abund_all,d["C+"]), "o", mfc="none", ms=ms,label="C+")
ax[2,1].plot(nH,getindex.(abund_all,d["C"]) , "o", mfc="none", ms=ms,label="C")
ax[2,1].plot(nH,getindex.(abund_all,d["CO"]), "o", mfc="none", ms=ms,label="CO")
ax[2,1].plot(nH,getindex.(abund_all,d["H+"]), "o", mfc="none", ms=ms,label="Hp")
ax[2,1].set_xscale("log")
ax[2,1].set_yscale("log")
ax[2,1].set_xlabel("nH [cm^-3]")
ax[2,1].set_ylabel("x_i")
ax[2,1].set_ylim(1e-9,2e-4)
ax[2,1].legend(loc="best", frameon=false, markerscale=5.0/ms)

ax[2,2].plot(NH,getindex.(abund_all,d["C+"]), "o", mfc="none", ms=ms,label="C+")
ax[2,2].plot(NH,getindex.(abund_all,d["C"]) , "o", mfc="none", ms=ms,label="C")
ax[2,2].plot(NH,getindex.(abund_all,d["CO"]), "o", mfc="none", ms=ms,label="CO")
ax[2,2].plot(NH,getindex.(abund_all,d["H+"]), "o", mfc="none", ms=ms,label="Hp")
ax[2,2].set_xscale("log")
ax[2,2].set_yscale("log")
ax[2,2].set_xlabel("NH_tot [cm^-2]")
ax[2,2].set_ylabel("x_i")
ax[2,2].set_ylim(1e-9,2e-4)
ax[2,2].legend(loc="best", frameon=false, markerscale=5.0/ms)

savefig(fname*".png")

#using Statistics
#const fac_col = 8.86780424e23
#med_col_tot = dropdims(median(col_tot, dims=1), dims=1) .* fac_col
#med_col_H2 = dropdims(median(col_H2, dims=1), dims=1) .* fac_col
#med_col_CO = dropdims(median(col_CO, dims=1), dims=1) .* fac_col
#xH2 = abund[1,:];
#xHp = abund[2,:];
#CO = abund[3,:];

#using Healpix
#const NSIDE=4
#m = Map{Float64, RingOrder}(NSIDE);
#m.pixels[:] = log10.(column_all .+ 1e-2*minimum(column_all[column_all.!=0.0]));
