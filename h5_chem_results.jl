using Serialization
using StaticArrays
using PyPlot
using HDF5
using Random
using Statistics
using StatsBase

ENV["MPLBACKEND"] = "Agg"
push!(LOAD_PATH, pwd())
using GadgetReader

#file_path = "/Users/chu/simulations/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Lsh100"

#file_path = "/ptmp/huchiayu/snapshots/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Lsh50"
#file_path = "/ptmp/huchiayu/snapshots/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Lsh100"
#const Zp = 1.0
file_path = "/ptmp/huchiayu/snapshots/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Z0p3"
const Zp = 0.3
#file_path = "/ptmp/huchiayu/snapshots/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Z0p1"
#const Zp = 0.1

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

const fac_col = (UnitMass_in_g/UnitLength_in_cm^2)/PROTONMASS
const facNHtoAv = 5.35e-22

const noneq="eq"
#const noneq="neq"

const T = Float64

function main()

    snaps = collect(150:10:1000)
    #snaps = collect(640:640)
    Ntime = length(snaps)
    cc=0 #counter for arrays
    Nbin = 40

    xH_n_time  = zeros(Nbin,Ntime)
    xH2_n_time = zeros(Nbin,Ntime)
    xCp_n_time = zeros(Nbin,Ntime)
    xC_n_time  = zeros(Nbin,Ntime)
    xCO_n_time = zeros(Nbin,Ntime)
    xH_N_time  = zeros(Nbin,Ntime)
    xH2_N_time = zeros(Nbin,Ntime)
    xCp_N_time = zeros(Nbin,Ntime)
    xC_N_time  = zeros(Nbin,Ntime)
    xCO_N_time = zeros(Nbin,Ntime)
    NH_n_time = zeros(Nbin,Ntime)
    NH2_n_time = zeros(Nbin,Ntime)
    NCO_NH2_time = zeros(Nbin,Ntime)


    for itime in snaps
        cc += 1
        snap = ""
        if itime < 10
            snap = "00" * string(itime)
        elseif itime < 100
            snap = "0" * string(itime)
        else
            snap = string(itime)
        end


Npart, pos, vel, rho, u, mass, hsml, id_gas,
    abund, fH2, fdust, col_tot, col_H2, col_CO, Tdust,
    N_star, pos_star, vel_star, m_star, sftime, id_star,
    boxsize, time_sim = read_snap(file_path * "/snap_" * snap * ".hdf5");


    println("=============== snapshot ", itime, " ===============")
    println("load data from file...")
if noneq=="neq"
    fname_base = file_path * "/chem-neqH2Hp-noCOic-TF3-" * snap * "-2"
elseif noneq=="eq"
    fname_base = file_path * "/chem-eqH2Hp-noCOic-TF3-" * snap * "-2"
end
fname = fname_base * ".hdf5"
header = h5readattr(fname, "/Header")
all_species::Vector{String} = header["all_species"]
time_chem::T   = header["time"]
facSFR::T      = header["facSFR"]
abund_all::Matrix{T}   = h5read(fname, "Chemistry/Abundances");
id_chem::Vector{Int}     = h5read(fname, "Chemistry/ID");
NH_eff::Vector{T}      = h5read(fname, "Chemistry/NH_eff");
NH2_eff::Vector{T}     = h5read(fname, "Chemistry/NH2_eff");
NCO_eff::Vector{T}    = h5read(fname, "Chemistry/NCO_eff");
println("done!")

#abund_all .+= 1e-15

N_spec = size(abund_all)[1]
d = Dict(all_species .=> collect(1:N_spec) );
d[""] = d["PHOTON"] = d["CRP"] = d["CRPHOT"] = 0


xH2 = abund_all[d["H2"], :]
xH  = abund_all[d["H"], :]
xHp = abund_all[d["H+"], :]
xCp = abund_all[d["C+"], :]
xC  = abund_all[d["C"], :]
xCO = abund_all[d["CO"], :]
xelec = abund_all[d["e-"], :]
xHe = abund_all[d["He"], :]

xH2[xH2.<0].=0
xH[xH.<0].=0
xHp[xHp.<0].=0
xCp[xCp.<0].=0
xC[xC.<0].=0
xCO[xCO.<0].=0
xHe[xHe.<0].=0
xelec[xelec.<0].=0

xHe
mu = @. 1 / (XH * (xH + xH2 + xHp + xHe + xelec))

nH = rho .* (XH * UnitDensity_in_cgs / PROTONMASS)
temp = u./mu ./ (1.5 * BOLTZMANN /PROTONMASS / 1e10)

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
        #@show y_median[i], y_low[i], y_high[i]
    end
    #ax.plot(x_center, y_mean, lsty, color=color, label=lbl)
    #ax.fill_between(x_center, y_mean.-y_std, y_mean.+y_std, color=color, alpha=alpha)
    ax.plot(x_center, y_median, lsty, color=color, label=lbl)
    ax.fill_between(x_center, y_low, y_high, color=color, alpha=alpha)
    return x_center, y_median
end

clf()
FS=17
#fig, ax = PyPlot.subplots(1, 2, figsize=(8.5,4))
fig, axes = PyPlot.subplots(2, 3, figsize=(14,8))
al = 0.3
nmin, nmax = 1, 5
Nmin,Nmax=20, 23
NH2min,NH2max=18,23

ax = axes[1,1]
nH_x, xH_n  = plot_xy_binned(log10.(nH), log10.(xH), Nbin, nmin, nmax, "blue", "-", ax, al, L"{\rm H}")
nH_x, xH2_n = plot_xy_binned(log10.(nH), log10.(2 .*xH2), Nbin, nmin, nmax, "red", "-", ax, al, L"{\rm H_2(2x)}")
ax.set_xlim(1, nmax)
ax.set_ylim(-2.0, 0.1)
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="lower right", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

#ymin,ymax=-7,-3.5
ymin,ymax=-7+log10(Zp),-3.7+log10(Zp)
ax = axes[2,1]
nH_x, xCp_n = plot_xy_binned(log10.(nH), log10.(xCp), Nbin, nmin, nmax, "blue" , "-", ax, al, L"{\rm C^+}")
nH_x, xC_n  = plot_xy_binned(log10.(nH), log10.(xC) , Nbin, nmin, nmax, "green", "-", ax, al, L"{\rm C}")
nH_x, xCO_n = plot_xy_binned(log10.(nH), log10.(xCO), Nbin, nmin, nmax, "red"  , "-", ax, al, L"{\rm CO}")
#plot_xy_binned(log10.(nH), log10.(abund[3,:].*(1.4/3)), Nbin, nmin, nmax, "red"  , "--", ax, al, L"{\rm CO}")
ax.set_xlim(1, nmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="upper left", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

ax = axes[1,2]
NH_x, xH_N  = plot_xy_binned(log10.(NH_eff), log10.(xH), Nbin, Nmin, Nmax, "blue", "-", ax, al, L"{\rm H}")
NH_x, xH2_N = plot_xy_binned(log10.(NH_eff), log10.(2 .*xH2), Nbin, Nmin, Nmax, "red", "-", ax, al, L"{\rm H_2(2x)}")
ax.set_xlim(Nmin,Nmax)
ax.set_ylim(-2.0, 0.1)
ax.set_xlabel(L"log_{10} N_H\ [{\rm cm^{-2}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="lower right", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

ax = axes[2,2]
NH_x, xCp_N = plot_xy_binned(log10.(NH_eff), log10.(xCp), Nbin, Nmin, Nmax, "blue" , "-", ax, al, L"{\rm C^+}")
NH_x, xC_N  = plot_xy_binned(log10.(NH_eff), log10.(xC) , Nbin, Nmin, Nmax, "green", "-", ax, al, L"{\rm C}")
NH_x, xCO_N = plot_xy_binned(log10.(NH_eff), log10.(xCO), Nbin, Nmin, Nmax, "red"  , "-", ax, al, L"{\rm CO}")
ax.set_xlim(Nmin,Nmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(L"log_{10} N_H\ [{\rm cm^{-2}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} x_i", fontsize=FS+4)
ax.legend(loc="upper left", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

#NH_sim = [-log(mean(exp.(-col_tot[:,i].*facNHtoAv.*fac_col.*XH))) / facNHtoAv for i in 1:Npart]
#NH2_sim = [-log(mean(exp.(-col_H2[:,i].*facNHtoAv.*fac_col./2))) / facNHtoAv for i in 1:Npart]
#NCO_sim = [-log(mean(exp.(-col_CO[:,i].*facNHtoAv.*fac_col./28))) / facNHtoAv for i in 1:Npart]

ax = axes[1,3]
#plot_xy_binned(log10.(nH), log10.(NH_sim) , Nbin, nmin, nmax, "blue" , "-", ax, al, L"{\rm NH_{sim}}")
#plot_xy_binned(log10.(nH), log10.(NH2_sim), Nbin, nmin, nmax, "red" , ":", ax, al, L"{\rm NH2_{sim}}")
nH_x, NH_n  = plot_xy_binned(log10.(nH), log10.(NH_eff) , Nbin, nmin, nmax, "red"  , "-", ax, al, L"{\rm N_{H,eff}}")
nH_x, NH2_n = plot_xy_binned(log10.(nH), log10.(NH2_eff), Nbin, nmin, nmax, "blue"  , ":", ax, al, L"{\rm N_{H_2,eff}}")
ax.set_xlim(nmin,nmax)
ax.set_ylim(NH2min,NH2max)
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} N\ [{\rm cm^{-2}}]", fontsize=FS+2)
ax.legend(loc="lower right", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")

ax = axes[2,3]
#plot_xy_binned(log10.(nH), log10.(NCO_sim), Nbin, nmin, nmax, "blue" , "--", ax, al, L"{\rm NCO_{sim}}")
NH2_x, NCO_NH2 = plot_xy_binned(log10.(NH2_eff), log10.(NCO_eff), Nbin, NH2min, NH2max, "blue"  , "-", ax, al, L"{\rm NCO_{eff}}")
ax.plot([NH2min,NH2max], [NH2min,NH2max].+log10(1.4e-4*Zp/0.5), ":", c="grey")
ax.set_xlim(NH2min,NH2max)
ax.set_ylim(12,20)
ax.set_xlabel(L"log_{10} N_{H_2}\ [{\rm cm^{-2}}]", fontsize=FS+2)
ax.set_ylabel(L"log_{10} N_{CO}\ [{\rm cm^{-2}}]", fontsize=FS+2)
#ax.legend(loc="best", frameon=false, fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true, ls=":")
fig.tight_layout()
#savefig(fname_base*".png")
savefig(file_path*"/transition_"*noneq*"_"*snap*".png")
close(fig)
#return nH, temp, mu, NH_eff, NH2_eff, NCO_eff, xH2, xCO
xH_n_time[:,cc] = xH_n
xH2_n_time[:,cc] = xH2_n
xCp_n_time[:,cc] = xCp_n
xC_n_time[:,cc] = xC_n
xCO_n_time[:,cc] = xCO_n
xH_N_time[:,cc] = xH_N
xH2_N_time[:,cc] = xH2_N
xCp_N_time[:,cc] = xCp_N
xC_N_time[:,cc] = xC_N
xCO_N_time[:,cc] = xCO_N
NH_n_time[:,cc] = NH_n
NH2_n_time[:,cc] = NH2_n
NCO_NH2_time[:,cc] = NCO_NH2
fname = file_path * "/data_transition_" * noneq * ".hdf5"
fid=h5open(fname,"w")
grp_part = g_create(fid,"gas");
h5write(fname, "gas/xH_n_time"   , xH_n_time)
h5write(fname, "gas/xH2_n_time"  , xH2_n_time)
h5write(fname, "gas/xCp_n_time"  , xCp_n_time)
h5write(fname, "gas/xC_n_time"   , xC_n_time)
h5write(fname, "gas/xCO_n_time"  , xCO_n_time)
h5write(fname, "gas/xH_N_time"   , xH_N_time)
h5write(fname, "gas/xH2_N_time"  , xH2_N_time)
h5write(fname, "gas/xCp_N_time"  , xCp_N_time)
h5write(fname, "gas/xC_N_time"   , xC_N_time)
h5write(fname, "gas/xCO_N_time"  , xCO_N_time)
h5write(fname, "gas/NH_n_time"   , NH_n_time)
h5write(fname, "gas/NH2_n_time"  , NH2_n_time)
h5write(fname, "gas/NCO_NH2_time", NCO_NH2_time)
close(fid)

end #for-loop
end #main

#nH, temp, mu, NH_eff, NH2_eff, NCO_eff, xH2, xCO =
main()
#savefig(file_path * "/nH_vs_x_bin_" * noneq *".png")

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
