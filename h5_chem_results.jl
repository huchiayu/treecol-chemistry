using Serialization
using StaticArrays
using PyPlot
using HDF5
using Random
using Statistics
import StatsBase

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

    N_star::Int64 = header["NumPart_ThisFile"][5]

    pos_star::Matrix{T} = h5read(filename, "PartType4/Coordinates");
    vel_star::Matrix{T} = h5read(filename, "PartType4/Velocities");
    m_star::Vector{T}   = h5read(filename, "PartType4/Masses");
    sftime::Vector{T}      = h5read(filename, "PartType4/StellarFormationTime");
    id_star::Vector{Int64} = h5read(filename, "PartType4/ParticleIDs");

    return N_gas, pos_gas, vel_gas, rho, u, m_gas, hsml, id_gas,
        abund, fH2ss, fH2dust, col_tot, col_H2, col_CO, Tdust,
        N_star, pos_star, vel_star, m_star, sftime, id_star,
        boxsize, time
end

#function read_chem(filename)
    #return NH_eff
#end

#file_path = "../../../simulations/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft4_SFLJ4_eff0p5_stoIMFfix_rngSF_convSF"
#num = "550"
file_path = "../../../simulations/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Z0p1"
num = "415"

Npart, pos, vel, rho, u, mass, hsml, id_gas,
    abund, fH2, fdust, col_tot, col_H2, col_CO, Tdust,
    N_star, pos_star, vel_star, m_star, sftime, id_star,
    boxsize, time_sim = read_snap(file_path * "/snap_" * num * ".hdf5");

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

nH = rho .* (XH * UnitDensity_in_cgs / PROTONMASS)
temp = u ./ (1.5 * BOLTZMANN /PROTONMASS / 1e10)

# load from file
#fname = file_path * "/chem3D-neqH2Hp-550-5"
fname = file_path * "/chem3D-neqH2Hp-" * num * "-1.hdf5"
header = h5readattr(fname, "/Header")
all_species = header["all_species"]
time_chem   = header["time"]
facSFR      = header["facSFR"]
abund_all   = h5read(fname, "Chemistry/Abundances");
id_chem     = h5read(fname, "Chemistry/ID");
NH_eff      = h5read(fname, "Chemistry/NH_eff");
NH2_eff     = h5read(fname, "Chemistry/NH2_eff");
NCO_eff     = h5read(fname, "Chemistry/NCO_eff");

N_spec = size(abund_all)[1]
d = Dict(all_species .=> collect(1:N_spec) );
d[""] = d["PHOTON"] = d["CRP"] = d["CRPHOT"] = 0


xH2 = abund_all[d["H2"], :]
xH  = abund_all[d["H"], :]
xHp = abund_all[d["H+"], :]
xCp = abund_all[d["C+"], :]
xC  = abund_all[d["C"], :]
xCO = abund_all[d["CO"], :]


function plot_xy_binned(x_data, y_data, nbin, x_min, x_max, color, lsty, ax, alpha, lbl)
    #y_mean = zeros(nbin) .- Inf
    #y_std  = zeros(nbin) .- Inf
    y_median = zeros(nbin) .- Inf
    y_low  = zeros(nbin) .- Inf
    y_high = zeros(nbin) .- Inf
    dx = (x_max - x_min) / nbin
    x_edge = x_min .+ (collect(0:nbin))  .* dx
    x_center = @. 0.5 * (x_edge[2:end] + x_edge[1:end-1])
    for i in 1:nbin
        idx = (abs.(x_data .- x_center[i]) .< (0.5 * dx)) .& (y_data .!= -Inf)
        if sum(idx) == 0 continue end
        #y_mean[i] = mean(y_data[idx])
        #y_std[i] = std(y_data[idx])
        y_median[i] = median(y_data[idx])
        y_low[i] = StatsBase.percentile(y_data[idx], 16)
        y_high[i] = StatsBase.percentile(y_data[idx], 84)
    end
    #return x_center, y_mean, y_std, y_median, y_low, y_high

    #ax.plot(x_center, y_mean, lsty, color=color, label=lbl)
    #ax.fill_between(x_center, y_mean.-y_std, y_mean.+y_std, color=color, alpha=alpha)
    ax.plot(x_center, y_median, lsty, color=color, label=lbl)
    ax.fill_between(x_center, y_low, y_high, color=color, alpha=alpha)
end

clf()
FS=17
#fig, ax = PyPlot.subplots(1, 2, figsize=(8.5,4))
fig, axes = PyPlot.subplots(2, 2, figsize=(10,8))
Nbin = 30
alpha = 0.3
nmin, nmax = 0.2, 5.
ax = axes[1,1]
plot_xy_binned(log10.(nH), log10.(xH), Nbin, nmin, nmax, "blue", "-", ax, alpha, L"{\rm H}")
plot_xy_binned(log10.(nH), log10.(2 .*xH2), Nbin, nmin, nmax, "red", "-", ax, alpha, L"{\rm H_2(2x)}")
ax.set_xlim(1, nmax)
ax.set_ylim(-2.0, 0.1)
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="best", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

#ymin,ymax=-7,-3.5
ymin,ymax=-8,-4.5
ax = axes[2,1]
plot_xy_binned(log10.(nH), log10.(xCp), Nbin, nmin, nmax, "blue" , "-", ax, alpha, L"{\rm C^+}")
plot_xy_binned(log10.(nH), log10.(xC) , Nbin, nmin, nmax, "green", "-", ax, alpha, L"{\rm C}")
plot_xy_binned(log10.(nH), log10.(xCO), Nbin, nmin, nmax, "red"  , "-", ax, alpha, L"{\rm CO}")
ax.set_xlim(1, nmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="upper left", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

Nmin,Nmax=20.8,23
ax = axes[1,2]
plot_xy_binned(log10.(NH_eff), log10.(xH), Nbin, Nmin, Nmax, "blue", "-", ax, alpha, L"{\rm H}")
plot_xy_binned(log10.(NH_eff), log10.(2 .*xH2), Nbin, Nmin, Nmax, "red", "-", ax, alpha, L"{\rm H_2(2x)}")
#ax.set_xlim(1, nmax)
ax.set_ylim(-2.0, 0.1)
ax.set_xlabel(L"log_{10} N_H\ [{\rm cm^{-2}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="best", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

ax = axes[2,2]
plot_xy_binned(log10.(NH_eff), log10.(xCp), Nbin, Nmin, Nmax, "blue" , "-", ax, alpha, L"{\rm C^+}")
plot_xy_binned(log10.(NH_eff), log10.(xC) , Nbin, Nmin, Nmax, "green", "-", ax, alpha, L"{\rm C}")
plot_xy_binned(log10.(NH_eff), log10.(xCO), Nbin, Nmin, Nmax, "red"  , "-", ax, alpha, L"{\rm CO}")
#ax.set_xlim(1, nmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(L"log_{10} N_H\ [{\rm cm^{-2}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="upper left", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")
fig.tight_layout()
savefig(fname*"_binned.png")

#=
nH = rho .* 287
temp = 80 .* u
idx = (nH.>1) .& (temp.<3e3)

nH = nH[idx]
NH = NH[idx]
abund_all = abund_all[idx]

#boxsize = 0.15
#r = [sqrt(sum((X[i].- 0.5*boxsize).^2)) for i in eachindex(X)]
=#
