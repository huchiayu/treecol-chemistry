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

#file_path = "/ptmp/huchiayu/snapshots/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Lsh500"
#file_path = "/ptmp/huchiayu/snapshots/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Lsh50"
#file_path = "/ptmp/huchiayu/snapshots/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Lsh100"
#const Zp = 1.0
#file_path = "/ptmp/huchiayu/snapshots/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Z3"
#const Zp = 3.0
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


#const noneq="neq"
const noneq="eq"

const T = Float64

function main()
snaps = collect(150:10:1000)
#snaps = collect(770:770)
Ntime = length(snaps)
Nbin=101

h_gas_time=zeros(Nbin-1,Nbin-1,Ntime);
pdfx_gas_time=zeros(Nbin-1,Ntime);
pdfy_gas_time=zeros(Nbin-1,Ntime);
pdfx_H2_time=zeros(Nbin-1,Ntime);
pdfy_H2_time=zeros(Nbin-1,Ntime);
pdfx_CO_time=zeros(Nbin-1,Ntime);
pdfy_CO_time=zeros(Nbin-1,Ntime);
tsim = zeros(Ntime);
M_H2_time = zeros(Ntime);
M_CO_time = zeros(Ntime);
SFR_time1 = zeros(Ntime);
SFR_time2 = zeros(Ntime);
cc=0 #counter for arrays

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


mu = @. 1 / (XH * (xH + xH2 + xHp + xHe + xelec))
nH = rho .* (XH * UnitDensity_in_cgs / PROTONMASS)
#T_over_mu = u ./ (1.5 * BOLTZMANN /PROTONMASS / 1e10)
temp = u.*mu .* (1e10*PROTONMASS/(1.5*BOLTZMANN))


xmin,xmax,ymin,ymax=-4,6,0.5,8
xbins, ybins = range(xmin,xmax,length=Nbin), range(ymin,ymax,length=Nbin);
h_gas = fit(Histogram, (log10.(nH), log10.(temp)), (xbins, ybins) );
h_H2 = fit(Histogram, (log10.(nH), log10.(temp)), weights(xH2.*2), (xbins, ybins) );
h_CO = fit(Histogram, (log10.(nH), log10.(temp)), weights(xCO.*28), (xbins, ybins) );
xbin = 0.5*(h_gas.edges[1][1:end-1] + h_gas.edges[1][2:end])
ybin = 0.5*(h_gas.edges[2][1:end-1] + h_gas.edges[2][2:end])
pdfx_gas = dropdims(sum(h_gas.weights,dims=2), dims=2)
pdfy_gas = dropdims(sum(h_gas.weights,dims=1), dims=1)
pdfx_H2 = dropdims(sum(h_H2.weights,dims=2), dims=2)
pdfy_H2 = dropdims(sum(h_H2.weights,dims=1), dims=1)
pdfx_CO = dropdims(sum(h_CO.weights,dims=2), dims=2)
pdfy_CO = dropdims(sum(h_CO.weights,dims=1), dims=1)


tsim[cc] = time_sim
M_H2_time[cc] = sum(xH2.*2) * XH * mean(mass) * 1e10
M_CO_time[cc] = sum(xCO.*28) * XH * mean(mass) * 1e10

sfr_dt1 = 1e-2
sfr_dt2 = 3e-2
idx = @. ( time_sim > sftime ) & ( time_sim - sftime < sfr_dt1 )
SFR_time1[cc] = 1e10 * sum(abs.(m_star[idx])) / (sfr_dt1*1e9)
idx = @. ( time_sim > sftime ) & ( time_sim - sftime < sfr_dt2 )
SFR_time2[cc] = 1e10 * sum(abs.(m_star[idx])) / (sfr_dt2*1e9)


h_gas_time[:,:, cc] = h_gas.weights
pdfx_gas_time[:, cc] = pdfx_gas
pdfy_gas_time[:, cc] = pdfy_gas
pdfx_H2_time[:, cc] = pdfx_H2
pdfy_H2_time[:, cc] = pdfy_H2
pdfx_CO_time[:, cc] = pdfx_CO
pdfy_CO_time[:, cc] = pdfy_CO

FS=14
lw=2
abs_C = 1.4e-4 * 28

clf()
fig, axx = PyPlot.subplots(2, 2, figsize=(13,10))
ax = axx[1,1]
ax.imshow(log10.(transpose(h_gas.weights)), origin="lower", extent=(xmin, xmax, ymin, ymax), aspect="equal")
#ax.set_xlim(xmin,xmax)
#ax.set_ylim(ymin,ymax)
ax.axis([xmin,xmax,ymin,ymax])
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS)
ax.set_ylabel(L"log_{10} T\ [{\rm K}]", fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true,linestyle=":")

ax = axx[1,2]
ax.plot(tsim, M_H2_time, "-" , lw=lw, label=L"{\rm M_{H_2}} [M_\odot]")
ax.plot(tsim, M_CO_time.*1e3, "-" , lw=lw, label=L"{\rm M_{CO}\times 10^3} [M_\odot]")
ax.plot(tsim, SFR_time1.*1e9, ":" , lw=lw, label=L"{\rm SFR_1} [M_\odot Gyr^{-1}]")
ax.plot(tsim, SFR_time2.*1e9, ":" , lw=lw, label=L"{\rm SFR_2} [M_\odot Gyr^{-1}]")

ax.set_yscale("log")
ax.legend(loc="upper left", frameon=false, fontsize=FS, ncol=2)
ax.set_xlabel(L"time [Gyr]", fontsize=FS)
ax.set_ylabel(L"{\rm mass [M_\odot]}", fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.set_xlim(0,1)
ax.set_ylim(1e3,3e6)

ax = axx[2,1]
ax.plot(xbin, pdfx_gas, "-" , lw=lw, label=L"{\rm total}")
ax.plot(xbin, pdfx_H2 , "--", lw=lw, label=L"{\rm H_2}")
ax.plot(xbin, pdfx_CO./abs_C , ":" , lw=lw, label=L"{\rm CO/xC}")
ax.set_yscale("log")
ax.set_ylim(1, 1e5)
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS)
ax.set_ylabel(L"{\rm mass-weighted\ PDF}", fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true,linestyle=":")
ax.legend(loc="lower left", frameon=false, fontsize=FS)

ax = axx[2,2]
ax.plot(ybin, pdfy_gas, "-" , lw=lw, label=L"{\rm total}")
ax.plot(ybin, pdfy_H2 , "--", lw=lw, label=L"{\rm H_2}")
ax.plot(ybin, pdfy_CO./abs_C , ":" , lw=lw, label=L"{\rm CO/xC}")
ax.set_yscale("log")
ax.set_ylim(1, 1e5)
ax.set_xlabel(L"log_{10} T\ [{\rm K}]", fontsize=FS)
ax.set_ylabel(L"{\rm mass-weighted\ PDF}", fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true,linestyle=":")
ax.legend(loc="upper right", frameon=false, fontsize=FS)

savefig(file_path * "/PD_" * noneq * "_" * snap * ".png")
    close(fig)

    fname = file_path * "/data_PD_" * noneq * ".hdf5"
    fid=h5open(fname,"w")
    grp_part = g_create(fid,"gas");
    h5write(fname, "gas/h_gas_time"     , h_gas_time)
    h5write(fname, "gas/pdfx_gas_time"  , pdfx_gas_time)
    h5write(fname, "gas/pdfy_gas_time"  , pdfy_gas_time)
    h5write(fname, "gas/pdfx_H2_time"  , pdfx_H2_time)
    h5write(fname, "gas/pdfy_H2_time"  , pdfy_H2_time)
    h5write(fname, "gas/pdfx_CO_time"  , pdfx_CO_time)
    h5write(fname, "gas/pdfy_CO_time"  , pdfy_CO_time)
    h5write(fname, "gas/tsim"      , tsim)
    h5write(fname, "gas/M_H2_time" , M_H2_time)
    h5write(fname, "gas/M_CO_time" , M_CO_time)
    h5write(fname, "gas/SFR_time1" , SFR_time1)
    h5write(fname, "gas/SFR_time2" , SFR_time2)
    close(fid)

end #for-loop
end #main

main()
#snaps = collect(150:10:1000)
#snaps = 770:770
#nH, temp, mu, abund_all, d, xH2, xCO, time_sim = main(770)
#for itime in snaps
#    main(itime)
#end
