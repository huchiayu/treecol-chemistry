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

    return N_gas, pos_gas, vel_gas, rho, u, m_gas, hsml, id_gas,
        abund, fH2ss, fH2dust, col_tot, col_H2, col_CO, Tdust, boxsize, time
end

#num = "1000"
#file_path = "../../../simulations/isocloud_N1e4/"
file_path = "../../../simulations/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft4_SFLJ4_eff0p5_stoIMFfix_rngSF_convSF"

#Ngas, pos, vel, rho, u, mass, hsml, id,
#abund, fH2ss, fH2dust, col_tot, col_H2, col_CO, Tdust, boxsize, time =
#        read_snap(file_path * "/snap_" * num * ".hdf5");

# load from file
fname = file_path * "/chem3D-neqH2Hp-550-5"

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
temp = 80 .* u
idx = (nH.>1) .& (temp.<3e3)

nH = nH[idx]
NH = NH[idx]
abund_all = abund_all[idx]

#boxsize = 0.15
#r = [sqrt(sum((X[i].- 0.5*boxsize).^2)) for i in eachindex(X)]

function plot_xy_binned(x_data, y_data, nbin, x_min, x_max, color, lsty, ax, alpha, lbl)
    y_mean = zeros(nbin)
    y_median = zeros(nbin)
    y_low  = zeros(nbin)
    y_high = zeros(nbin)
    y_std  = zeros(nbin)
    dx = (x_max - x_min) / nbin
    x_edge = x_min .+ (collect(0:nbin))  .* dx
    x_center = @. 0.5 * (x_edge[2:end] + x_edge[1:end-1])
    for i in 1:nbin
        idx = abs.(x_data .- x_center[i]) .< (0.5 * dx)
        y_mean[i] = mean(y_data[idx])
        y_median[i] = median(y_data[idx])
        y_std[i] = std(y_data[idx])
        y_low[i] = StatsBase.percentile(y_data[idx], 16)
        y_high[i] = StatsBase.percentile(y_data[idx], 84)
    end
    #return x_center, y_mean, y_std, y_median, y_low, y_high

    #plot(x_center, y_mean, lsty, color=color, label=lbl)
    #fill_between(x_center, y_mean.-y_std, y_mean.+y_std, color=color, alpha=alpha)
    ax.plot(x_center, y_median, lsty, color=color, label=lbl)
    ax.fill_between(x_center, y_low, y_high, color=color, alpha=alpha)
end

ms = 0.03

clf()
FS=17
#fig, ax = PyPlot.subplots(1, 2, figsize=(8.5,4))
fig, axes = PyPlot.subplots(2, 2, figsize=(10,8))
Nbin = 30
alpha = 0.3
nmin, nmax = 0, 4.5
ax = axes[1,1]
plot_xy_binned(log10.(nH), log10.(getindex.(abund_all,d["H"])), Nbin, nmin, nmax, "blue", "-", ax, alpha, L"{\rm H}")
plot_xy_binned(log10.(nH), log10.(2 .*getindex.(abund_all,d["H2"])), Nbin, nmin, nmax, "red", "-", ax, alpha, L"{\rm H_2(2x)}")
ax.set_xlim(1, nmax)
ax.set_ylim(-2.0, 0.1)
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="best", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

ax = axes[2,1]
plot_xy_binned(log10.(nH), log10.(getindex.(abund_all,d["C+"])), Nbin, nmin, nmax, "blue" , "-", ax, alpha, L"{\rm C^+}")
plot_xy_binned(log10.(nH), log10.(getindex.(abund_all,d["C"])) , Nbin, nmin, nmax, "green", "-", ax, alpha, L"{\rm C}")
plot_xy_binned(log10.(nH), log10.(getindex.(abund_all,d["CO"])), Nbin, nmin, nmax, "red"  , "-", ax, alpha, L"{\rm CO}")
ax.set_xlim(1, nmax)
ax.set_ylim(-7,-3.9)
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="upper left", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

Nmin,Nmax=20.8,22.5
ax = axes[1,2]
plot_xy_binned(log10.(NH), log10.(getindex.(abund_all,d["H"])), Nbin, Nmin, Nmax, "blue", "-", ax, alpha, L"{\rm H}")
plot_xy_binned(log10.(NH), log10.(2 .*getindex.(abund_all,d["H2"])), Nbin, Nmin, Nmax, "red", "-", ax, alpha, L"{\rm H_2(2x)}")
#ax.set_xlim(1, nmax)
ax.set_ylim(-2.0, 0.1)
ax.set_xlabel(L"log_{10} N_H\ [{\rm cm^{-2}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="best", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

ax = axes[2,2]
plot_xy_binned(log10.(NH), log10.(getindex.(abund_all,d["C+"])), Nbin, Nmin, Nmax, "blue" , "-", ax, alpha, L"{\rm C^+}")
plot_xy_binned(log10.(NH), log10.(getindex.(abund_all,d["C"])) , Nbin, Nmin, Nmax, "green", "-", ax, alpha, L"{\rm C}")
plot_xy_binned(log10.(NH), log10.(getindex.(abund_all,d["CO"])), Nbin, Nmin, Nmax, "red"  , "-", ax, alpha, L"{\rm CO}")
#ax.set_xlim(1, nmax)
ax.set_ylim(-7,-3.9)
ax.set_xlabel(L"log_{10} N_H\ [{\rm cm^{-2}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="upper left", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")
fig.tight_layout()
savefig(fname*"_binned.png")
#=
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
#ax[2,1].plot(nH,getindex.(abund_all,d["H+"]), "o", mfc="none", ms=ms,label="Hp")
ax[2,1].set_xscale("log")
ax[2,1].set_yscale("log")
ax[2,1].set_xlabel("nH [cm^-3]")
ax[2,1].set_ylabel("x_i")
ax[2,1].set_ylim(1e-9,2e-4)
ax[2,1].legend(loc="best", frameon=false, markerscale=5.0/ms)

ax[2,2].plot(NH,getindex.(abund_all,d["C+"]), "o", mfc="none", ms=ms,label="C+")
ax[2,2].plot(NH,getindex.(abund_all,d["C"]) , "o", mfc="none", ms=ms,label="C")
ax[2,2].plot(NH,getindex.(abund_all,d["CO"]), "o", mfc="none", ms=ms,label="CO")
#ax[2,2].plot(NH,getindex.(abund_all,d["H+"]), "o", mfc="none", ms=ms,label="Hp")
ax[2,2].set_xscale("log")
ax[2,2].set_yscale("log")
ax[2,2].set_xlabel("NH_tot [cm^-2]")
ax[2,2].set_ylabel("x_i")
ax[2,2].set_ylim(1e-9,2e-4)
ax[2,2].legend(loc="best", frameon=false, markerscale=5.0/ms)

savefig(fname*".png")
=#

#const fac_col = 8.86780424e23
#med_col_tot = dropdims(median(col_tot, dims=1), dims=1) .* fac_col
#med_col_H2 = dropdims(median(col_H2, dims=1), dims=1) .* fac_col
#med_col_CO = dropdims(median(col_CO, dims=1), dims=1) .* fac_col
#xH2 = abund[1,:];
#xHp = abund[2,:];
#CO = abund[3,:];
