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
using PyPlot

#import Plots

#include("OctTree.jl")
#using OctTree
#include("ChemistryNetwork.jl")
using ChemistryNetwork



T = Float64
const N = 3


function test()
    nH = 100.
    temp = 50.0 #check!
    ξ = 1.3e-17 #H2
    IUV = 1.0
    Zp = 1.0
    abund = zeros(N_spec)
    dtime = 1e8
    @time sol, rrates = solve_equilibrium_abundances(abund, dtime, nH, temp, 2.04e21, 8.75e19, 1.37e15, 8.28e16, ξ, IUV, Zp)
end
x,rrates=test();

const Nbin = 200
const xmin = 16.
const xmax = 23.


xl = "NH [cm^-2]"


function runPDR()
    nH = 100.
    temp = 50.0 #check!
    ξ = 1.3e-17 #H2
    IUV = 1.0
    Zp = 1.0
    dtime = 1e8

    xmin = 16.
    xmax = 23.

    ab_vs_x = Vector{Vector{Float64}}(undef, Nbin)
    for i in 1:Nbin
        ab_vs_x[i] = zeros(N_spec)
    end
    d_lnNH = (xmax-xmin) / (Nbin-1)
    xbin = @. 10^(xmin + collect(0:Nbin-1) * d_lnNH)
    NHbin = xbin
    dNH = NHbin[2:Nbin] .- NHbin[1:Nbin-1]

    NH2bin = zeros(Nbin)
    NCObin = zeros(Nbin)
    NCbin = zeros(Nbin)

    reaction_rates=zeros(Nbin, N_reac)

    @time for i in eachindex(NHbin)
        xH2 = getindex.(ab_vs_x, iH2)
        xCO = getindex.(ab_vs_x, iCO)
        xC  = getindex.(ab_vs_x, iC )
        NH = NHbin[i]
        NH2 = NCO = NC = 1e10
        if i > 1
            NH2 = sum(xH2[1:i-1] .* dNH[1:i-1])
            NCO = sum(xCO[1:i-1] .* dNH[1:i-1])
            NC  = sum(xC[1:i-1]  .* dNH[1:i-1])
            #NH2 = sum(xH2[1:i-1] .* NHbin[1:i-1]) * d_lnNH
            #NCO = sum(xCO[1:i-1] .* NHbin[1:i-1]) * d_lnNH
        end
        NH2bin[i] = NH2
        NCObin[i] = NCO
        NCbin[i] = NC
        #@show i, NH, NH2,NCO, NC
        print(i," ")


        G = 2.8e-5 # for Zp = 1
        αG = IUV * kH2diss / kdust / nH * G
        i==1 ? println("αG/2 = ", αG/2) : nothing

        #xneq[1] = xH2_eq[i]
        #xneq[2] = xHp_eq[i]

        #abund_eq, reaction_rates[i,:] =
        solve_equilibrium_abundances(ab_vs_x[i], dtime, nH, temp, NH, NH2, NCO, NC, ξ, IUV, Zp)

        #abund_eq, reaction_rates[i,:] = solve_equilibrium_abundances(100., 100., 1e20, 1e18, 1e10, 1e12, ξ, IUV, Zp)

        calc_abund_derived(ab_vs_x[i], Zp)
        sumH  = sum( ab_vs_x[i] .* fac_H )
        sumC  = sum( ab_vs_x[i] .* fac_C )
        sumO  = sum( ab_vs_x[i] .* fac_O )
        #sumSi = sum( ab_vs_x[i] .* fac_Si )
        #sumS  = sum( ab_vs_x[i] .* fac_S )
        sumelec = sum( ab_vs_x[i] .* charge )
        (1.0 ≈ sumH)  ? nothing : error("sumH = " , sumH)
        (abC_s  * Zp ≈ sumC)  ? nothing : error("sumC = " , sumC)
        (abO_s  * Zp ≈ sumO)  ? nothing : error("sumO = " , sumO)
        #(abSi_s * Zp ≈ sumSi) ? nothing : error("sumSi = ", sumSi)
        #(abS_s  * Zp ≈ sumS)  ? nothing : error("sumS = " , sumS)
        (1.0 ≈ 1.0 + sumelec) ? nothing : error("sumelec = ",sumelec)

        #ab_vs_x[i] = abund_eq
    end

    xbin = @. 10^(xmin + (xmax-xmin) * collect(0:Nbin-1) / (Nbin-1));
    N_H = NHbin
    clf()
    #fig = figure("", figsize=(12,8))
    fig = figure("", figsize=(12,4))

    abC_tot  = [sum(ab_vs_x[i] .* fac_C)  for i in 1:Nbin]
    abO_tot  = [sum(ab_vs_x[i] .* fac_O)  for i in 1:Nbin]
    #abSi_tot = [sum(ab_vs_x[i] .* fac_Si) for i in 1:Nbin]
    #abS_tot  = [sum(ab_vs_x[i] .* fac_S)  for i in 1:Nbin]

    xminp, xmaxp = 16., 23.
    ax1 = subplot(131)
    plot(xbin, 2 .*getindex.(ab_vs_x, dict["H2"]), "-"  , label="H2")  #H2
    plot(xbin, getindex.(ab_vs_x, dict["H"]), "--" , label="H")  #H
    plot(xbin, getindex.(ab_vs_x, dict["e-"]), ":"  , label="e-")  #e-
    #plot(xbin, N_H, ":"  , label="e-")  #e-

    #plot(xbin, N_H2./N_H, "-."  , label="NH2/NH")
    #plot(xbin, N_CO./N_H, ":"  , label="NCO/NH")
    legend(loc="best", fontsize=9, ncol=1, frameon=false)
    xscale("log")
    yscale("log")
    #axis([10^xmin, 10^xmax, 1e-3, 2])
    axis([10.0^xminp, 10.0^xmaxp, 1e-6, 2])
    #xlabel("Z'")
    xlabel(xl)
    ylabel("x_i")
    grid(linestyle=":", linewidth=1)
    ax2 = ax1.twiny()
    Av = N_H .* 5.35e-22
    ax2.set_xlabel("Av")
    ax2.set_xlim((10.0^xminp, 10.0^xmaxp).*5.35e-22)
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    xminp, xmaxp = 20., 23.
    ax1 = subplot(132)
    #plot(Zbin, getindex.(ab_vs_Z, 5), "-")  #C+
    plot(xbin, getindex.(ab_vs_x, dict["H+"]), "-" , label="H+")  #H
    plot(xbin, getindex.(ab_vs_x, dict["H3+"]) , "--" , label="H3+")
    plot(xbin, getindex.(ab_vs_x, dict["He+"]) , "-." , label="He+")
    plot(xbin, getindex.(ab_vs_x, dict["O"])  , ":"  , label="O")   #O
    plot(xbin, getindex.(ab_vs_x, dict["OH"]) , "-" , label="OH")  #OH
    plot(xbin, getindex.(ab_vs_x, dict["O2"]) , "--" , label="O2")  #O2
    plot(xbin, getindex.(ab_vs_x, dict["H2O"]), "-." , label="H2O") #H2O
    legend(loc="best", fontsize=9, ncol=2, frameon=false)
    xscale("log")
    yscale("log")
    axis([10.0^xminp, 10.0^xmaxp, 1e-11, 1e-3])
    xlabel(xl)
    #ylabel("x_i")
    grid(linestyle=":", linewidth=1)
    ax2 = ax1.twiny()
    Av = N_H .* 5.35e-22
    ax2.set_xlabel("Av")
    ax2.set_xlim((10.0^xminp, 10.0^xmaxp).*5.35e-22)
    ax2.set_xscale("log")
    ax2.set_yscale("log")


    ax1 = subplot(133)
    #plot(Zbin, getindex.(ab_vs_Z, 5), "-")  #C+
    plot(xbin, getindex.(ab_vs_x, dict["C"])  , "-"  , label="C")   #C
    plot(xbin, getindex.(ab_vs_x, dict["CO"]) , "--" , label="CO")  #CO
    plot(xbin, getindex.(ab_vs_x, dict["C+"]) , "-." , label="C+")  #C+
    plot(xbin, getindex.(ab_vs_x, dict["CH"]) , ":" , label="CH")  #CH
    plot(xbin, getindex.(ab_vs_x, dict["CH2"]) , ":" , label="CH2")  #CH
    legend(loc="best", fontsize=9, ncol=2, frameon=false)
    xscale("log")
    yscale("log")
    axis([10.0^xminp, 10.0^xmaxp, 1e-11, 1e-3])
    xlabel(xl)
    grid(linestyle=":", linewidth=1)
    #ylabel("x_i")
    ax2 = ax1.twiny()
    Av = N_H .* 5.35e-22
    ax2.set_xlabel("Av")
    ax2.set_xlim((10.0^xminp, 10.0^xmaxp).*5.35e-22)
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    return ab_vs_x, NHbin, NH2bin, NCObin, NCbin, reaction_rates
end
ab_vs_x, N_H, N_H2, N_CO, N_C, rr = runPDR();
0






#=
function solve_chem_all_particles()

    #nH = 1000.
    #temp = 50.
    ξ = 1.3e-16 #H2
    IUV = 1.0
    Zp = 1.0
=#
