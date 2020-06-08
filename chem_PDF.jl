using Serialization
using StaticArrays
using PyPlot
using HDF5
using Random

push!(LOAD_PATH, pwd())
using GadgetReader


file_path = "../../../simulations/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft4_SFLJ4_eff0p5_stoIMFfix_rngSF_convSF"
num = "550"
Zp = 1.0
#file_path = "../../../simulations/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Z0p3"
#num = "475"
#Zp = 0.3
#file_path = "../../../simulations/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Z0p1"
#num = "415"
#Zp = 0.1
println("Zp = ", Zp)

Npart, pos, vel, rho, u, mass, hsml, id_gas,
    abund, fH2, fdust, col_tot, col_H2, col_CO, Tdust,
    N_star, pos_star, vel_star, m_star, sftime, id_star,
    boxsize, time_sim = read_snap(file_path * "/snap_" * num * ".hdf5");

# load from file
fname = file_path * "/chem3D-eqH2Hp-" * num * "-4.hdf5"
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
abund_all .+= 1e-15

nH = rho .* 287
temp = 80 .* u

fH_e  = abund_all[d["H"],:]
fHp_e = abund_all[d["H+"],:]
fH2_e = abund_all[d["H2"],:] .* 2;
fCO_e = abund_all[d["CO"],:] .* 28;
println("steady state H2:")
println("fCO/fH2 = ", sum(fCO_e[fCO_e.<2.9e-3]) / sum(fH2_e) )
println("M_H2 = ", sum(fH2_e) .* 10)
println("M_CO = ", sum(fCO_e) .* 10)

# load from file
fname = file_path * "/chem3D-neqH2Hp-" * num * "-4.hdf5"
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
abund_all .+= 1e-15
nH = rho .* 287
temp = 80 .* u

fH  = abund_all[d["H"],:]
fHp = abund_all[d["H+"],:]
fH2 = abund_all[d["H2"],:] .* 2;
fCO = abund_all[d["CO"],:] .* 28;
println("time-dependent H2:")
println("fCO/fH2 = ", sum(fCO[fCO.<2.9e-3]) / sum(fH2) )
println("M_H2 = ", sum(fH2) .* 10)
println("M_CO = ", sum(fCO) .* 10)


idx=@.( (fH2_e>0) & (nH>1.0) & (temp<3e3) );
#cc,xb,yb,im = hist2D(log10.(fH2[idx]), log10.(fH2[idx]./fH2_e[idx]), bins=[100,100], range=[[0, 4],[-5, 5]])
#imshow(transpose(log10.(cc)), origin="lower", extent=(0, 4, -5,5), aspect="auto")


nmin,nmax = -3,6
hGas,xx=hist(log10.(nH), bins=50, range=[nmin,nmax], histtype="step");
hH2eq,xx=hist(log10.(nH), bins=50, range=[nmin,nmax], histtype="step", weights=fH2_e);
hCOeq,xx=hist(log10.(nH), bins=50, range=[nmin,nmax], histtype="step", weights=fCO_e);
hH2neq,xx=hist(log10.(nH), bins=50, range=[nmin,nmax], histtype="step", weights=fH2);
hCOneq,xx=hist(log10.(nH), bins=50, range=[nmin,nmax], histtype="step", weights=fCO);
xx=(xx[2:end].+xx[1:end-1])./2;
clf()


#clf()
#fig, axes = PyPlot.subplots(1, 1, figsize=(5,4))
fig, ax = PyPlot.subplots(1, 1, figsize=(6,5))
lw=2
FS=14

#ax=axes[1]
#ax.set_title(L"{\rm time-dependent\ H_2}", fontsize=FS)
ax.set_title(L"Z = " * string(Zp) * L"Z_\odot", fontsize=FS)
ax.plot(xx, hGas  , "-", lw=lw, label=L"{\rm total\ gas}")
ax.plot(xx, hH2neq, "-", lw=lw, c="green", label=L"{\rm H_2}")
ax.plot(xx, hCOneq, "-", lw=lw, c="red", label=L"{\rm CO}")
ax.plot(xx, hH2eq , ":", lw=lw, c="green")
ax.plot(xx, hCOeq , ":", lw=lw, c="red")

ax.set_yscale("log")
ax.set_xlim(nmin,nmax)
ax.set_ylim(1e-3, 1e5)
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS)
ax.set_ylabel(L"{\rm mass-weighted\ PDF}", fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.grid(true,linestyle=":")
ax.legend(loc="best", frameon=false, fontsize=FS)

#ax=axes[2]
#ax.set_title(L"{\rm steady-state\ H_2}", fontsize=FS)
#ax.plot(xx, hGas  , "-", lw=lw, label=L"{\rm total\ gas}")
#=
ax.set_yscale("log")
ax.set_xlim(nmin,nmax)
ax.set_ylim(1e-3, 1e5)
ax.set_xlabel(L"log_{10} n_H\ [{\rm cm^{-3}}]", fontsize=FS)
ax.set_ylabel(L"{\rm mass-weighted\ PDF}", fontsize=FS)
ax.tick_params(labelsize=FS, axis="both")
ax.legend(loc="best", frameon=false, fontsize=FS)
ax.grid(true,linestyle=":")
=#
fig.tight_layout()

#savefig(file_path*"/PDF_H2_CO_"*num*".png")
