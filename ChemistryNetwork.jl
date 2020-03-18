module ChemistryNetwork

using DelimitedFiles
using Printf
using StaticArrays
using DifferentialEquations
#using OrdinaryDiffEq
#using NLsolve
#using PyPlot  #can't use Interact with it
#using Plots
#plotlyjs()
using Interpolations

using Sundials #for CVODE_BDF

export solve_equilibrium_abundances, calc_abund_derived, init_abund
export kH2diss, kdust
export iH2,iCO,iC,fac_H,fac_C,fac_O,dict

ssCO = readdlm("CO_shielding_functions/shield.03.5.69-557-36.dat"); #CO linewidth = 0.3, T_ex = 5
const Nbin_nco = 47
const Nbin_nh2 = 42
const Nco = SVector{Nbin_nco,Float64}(ssCO[9:55])
const Nh2 = SVector{Nbin_nh2,Float64}(ssCO[57:98])
const istart = 100 #start line number for 12CO

function get_fco(is)
    a=[]
    for i in 0:4
        a = vcat(a, ssCO[is+i,:])
    end
    a = a[1:Nbin_nco]
end

ssCOtable = zeros(Nbin_nco, Nbin_nh2)
for j in 1:Nbin_nh2
    jj = istart+(j-1)*5
    #@show j, jj
    ssCOtable[:, j] .= get_fco(jj)
end


fCOss = interpolate((Nco, Nh2), ssCOtable, Gridded(Linear()));
function fCOselfshield(Nco,Nh2)
    Nh2 = Nh2 >= 1e23 ? 1e23 : Nh2
    Nco = Nco >= 1e19 ? 1e19 : Nco
    Nh2 = Nh2 <= 1e10 ? 1e10 : Nh2
    Nco = Nco <= 1e10 ? 1e10 : Nco
    return fCOss(Nco, Nh2)
end

function fH2selfshield(Nh2)
    x = Nh2 * 2e-15
    b5 = 2.0
    return 0.965 / (1.0 + x/b5)^2 + 0.035 / sqrt(1 + x) * exp(-8.5e-4 * sqrt(1 + x))
end

function fCselfshield(NC,NH2)
    rH2 = 2.8e-22 * NH2
    #return 1.0
    #return exp( -1.6e-17 * NC )
    #return exp(-rH2) / (1.0 + rH2)
    return exp( -1.6e-17 * NC ) * exp(-rH2) / (1.0 + rH2)
end

function grain_recomb_H(Temp, psi)
    @assert (Temp > 0.0 && psi > 0.0)
    return 12.25e-14 / (1.0 + 8.074e-6 * psi^1.378 * (1.0 + 508.7*Temp^0.01586*psi^(-0.4723-1.102e-5*log(Temp))))
end

function grain_recomb_He(Temp, psi)
    return 5.572e-14 / (1.0 + 3.185e-7 * psi^1.512 * (1.0 + 5115*Temp^3.903e-7*psi^(-0.4956-5.494e-7*log(Temp))))
end

function grain_recomb_C(Temp, psi)
    return 45.58e-14 / (1.0 + 6.089e-3 * psi^1.128 * (1.0 + 433.1*Temp^0.04845*psi^(-0.8120-1.333e-4*log(Temp))))
end

function grain_recomb_Si(Temp, psi)
    return 2.166e-14 / (1.0 + 5.678e-8 * psi^1.874 * (1.0 + 43750*Temp^1.635e-6*psi^(-0.8964-7.538e-5*log(Temp))))
end

function read_umist_line(umist_line)
    num = umist_line[1]
    typ = umist_line[2]
    R::Vector{String} = [umist_line[3], umist_line[4]]
    P::Vector{String} = [umist_line[5], umist_line[6], umist_line[7], umist_line[8]]
    #alpha, beta, gamma, Tl, Tu, ST, ACC =
    #    umist_line[10], umist_line[11], umist_line[12],
    #    umist_line[13], umist_line[14], umist_line[15], umist_line[16]
    alpha, beta, gamma =
        umist_line[10], umist_line[11], umist_line[12]
    return num, typ, R, P, alpha, beta, gamma
end

function print_reaction(num, typ, R, P, alpha, beta, gamma)
    @printf "(%4s)" num
    #@printf "[%3s]" typ

    @printf " [α,β,γ = "
    @printf "%10.2e," alpha
    @printf "%10.2e," beta
    @printf "%10.2e]" gamma

    @printf "%8s" R[1]
    @printf "     + %8s" R[2]
    if typ == "GR"
        @printf "(gr)-->%8s" P[1]
    else
        @printf "    -->%8s" P[1]
    end
    if(P[2]!="")
        @printf "     + %8s" P[2]
    end
    if(P[3]!="")
        @printf "     + %8s" P[3]
    end
    if(P[4]!="")
        @printf "     + %8s" P[4]
    end
    print("\n")
end

function print_reaction(num, typ, R, P)
    @printf "%7s" R[1]
    @printf "     + %7s" R[2]
    if typ == "GR"
        @printf " (gr)-->%7s" P[1]
    else
        @printf "     -->%7s" P[1]
    end
    if(P[2]!="")
        @printf "     + %7s" P[2]
    end
    if(P[3]!="")
        @printf "     + %7s" P[3]
    end
    if(P[4]!="")
        @printf "     + %7s" P[4]
    end
    print("\n")
end

function search_species_reactants(network, s1::String)
    reac=[]
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        if any(s1 .== R)
            print_reaction(n, t, R, P, a, b, c)
            push!(reac, n)
        end
    end
    reac
end

function search_species_reactants(network, s1::String, s2::String)
    reac=[]
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        if any(s1 .== R) && any(s2 .== R)
            @printf "(i=%.4d) " i
            print_reaction(n, t, R, P, a, b, c)
            push!(reac, n)
        end
    end
    reac
end

function search_species_products(network, s1::String)
    reac=[]
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        if any(s1 .== P)
            print_reaction(n, t, R, P, a, b, c)
            push!(reac, n)
        end
    end
    reac
end

function search_species_products(network, s1::String, s2::String)
    reac=[]
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        if any(s1 .== P) && any(s2 .== P)
            print_reaction(n, t, R, P, a, b, c)
            push!(reac, n)
        end
    end
    reac
end

function print_all_reactions(network)
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        #@printf "(%3d)  " i
        print_reaction(n, t, R, P, a, b, c)
    end
end

function get_all_species_in_network(network)
    species::Vector{String} = []
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        #@printf "(%3d)  " i
        #print_reaction(n, t, R, P, a, b, c)
        [push!(species, R[j]) for j in eachindex(R)]
        [push!(species, P[k]) for k in eachindex(P)]
    end
    species = unique(species)
    filter!(s -> (s != "")&&(s != "PHOTON")&&(s != "CRPHOT")&&(s != "CRP"), species)
end

s="HSi3Mg2+"
function breakdown_species(s::String)
    charge::Int32 = 0
    elements = String[]
    counts   = Int32[]
    mult::Int32 = 1
    skipflag::Bool = false

    if s == "e-"
        return -1, elements, counts
    end

    for i in length(s):-1:1
        if skipflag
            #@show skipflag
            skipflag = false
            continue
        end
        char = s[i]
        #@show char
        isuppercase(char)
        if char == '+'
            charge += 1
        elseif char == '-'
            charge -= 1
        elseif isletter(char)
            elem = ""
            if islowercase(char)
                elem = s[i-1:i]
                skipflag = true
            else
                elem = string(char)
            end
            push!(elements, elem)
            push!(counts, mult)
            mult = 1
        else
            mult = parse(Int32, char)
        end
    end
    return charge, elements, counts
end
breakdown_species(s)


reactions = Int32[]

rea_H = [6155, 75, 4921, 6093, 489, 731, 730, 2604, 1278, 1277, 2895, 733, 406,
408, 434, 6162, 492]
rea_He = [734, 3431, 491, 458, 2677, 2610, 2678, 3065, 6094,6166, 3517, 6173, 967]
rea_O = [2686, 2951, 2689, 2662, 1296, 1295, 1293, 1294, 1266, 1264, 1265, 1427, 5643, 5199,
1646, 4029, 408, 1613, 2789, 2767, 2904, 378, 5196, 1643, 1642, 5705, 948, 405, 1421, 5412]
rea_CH = [2852, 2647, 6095, 2648, 1175, 1176, 1177, 1158, 1161, 1160, 1159, 3055, 5370, 5374]
rea_CO = [5199, 1646, 488, 2654, 1349, 5310, 5231, 5232, 64, 1349, 823, 6158, 707, 194, 876]
#rea_photo = [5997, 5915, 5989, 5887, 5888, 5902]
rea_suprathermal = [2618, 2189]

reactions = vcat(rea_H, rea_He, rea_O, rea_CH, rea_CO, rea_suprathermal)
unique!(reactions);

umist = readdlm("RATE12.dist.txt",':')
reduced_network = umist[reactions,:]

#const N_reac = size(reduced_network, 1)
const all_species = get_all_species_in_network(reduced_network);
#const N_spec = size(all_species,1);


#add photoprocesses automatically
charge_a = @. Int32(getindex(breakdown_species(all_species) , 1));
photo_reactions = []
for i in 1:size(all_species,1)
    #if charge_a[i] <= 0
    if true
        r = search_species_reactants(umist, "PHOTON", all_species[i])
        if all_species[i] != "CO" #will do CO shielding manually
            global photo_reactions = vcat(r, photo_reactions)
        end
    end
end
all_reactions = vcat(reactions, photo_reactions)
unique!(all_reactions);
reduced_network = umist[all_reactions,:];

const N_reac = size(reduced_network, 1)
const all_species = get_all_species_in_network(reduced_network);
println(all_species)
const N_spec = size(all_species,1)
@show N_reac, N_spec


const N_reac = size(reduced_network, 1)
const all_species = get_all_species_in_network(reduced_network);
print_all_reactions(reduced_network)
println(all_species)
const N_spec = size(all_species,1)
@show N_reac, N_spec

all_species_2 = vcat(all_species , ["", "PHOTON", "CRP", "CRPHOT"]);

reac_num=[]
for i in 1:size(umist,1)
    #println(i)
    n, t, R, P, a, b, c = read_umist_line(umist[i,:])
    if any(R[1].==all_species_2) && any(R[2].==all_species_2) && any(P[1].==all_species_2) && any(P[2].==all_species_2) && any(P[3].==all_species_2) && any(P[4].==all_species_2)
        #alpha_a[i] = 4.67e-12
        #gamma_a[i] = 3.76
        #@show a, b, c
        print_reaction(n, t, R, P, a, b, c)
        push!(reac_num,n)
    end
end

reduced_network = umist[reac_num,:]

#add 7 more reactions (H2 & CO dissociation & grain-processes)
more_reactions = readdlm("user-rates.txt",':')
reduced_network = vcat(reduced_network, more_reactions);
const N_reac = size(reduced_network, 1)


dict = Dict(all_species .=> collect(1:N_spec) );
dict[""] = dict["PHOTON"] = dict["CRP"] = dict["CRPHOT"] = 0

charge_a = @. Int32(getindex(breakdown_species(all_species) , 1))
elements = @. getindex(breakdown_species(all_species) , 2)
counts = @. getindex(breakdown_species(all_species) , 3)


unique_elements_a = String[]
for i in 1:N_spec
    for j in eachindex(elements[i])
        push!(unique_elements_a, elements[i][j])
    end
end
unique!(unique_elements_a)
const N_uni_elem = length(unique_elements_a)
const uni_elem = SVector{N_uni_elem}(unique_elements_a)
#N_conserve = length(unique_elements) + 1  #+1 for electron

#species_derived = ["e-", "H", "He", "C", "O", "Si", "S"]
species_derived = unique_elements_a
push!(species_derived, "e-")
@show species_derived


species_noneq = []
#species_noneq = ["H2"]
#xneq = [0.0]
#species_noneq = ["H2", "H+"]
#xneq = [0.0, 0.0]

const N_neq = length(species_noneq)
const NONEQ = N_neq > 0

if NONEQ
    const ineq = [dict[species_noneq[i]] for i in 1:N_neq]

    for i in 1:N_neq
        push!(species_derived, species_noneq[i])
    end
end
#@show species_derived
#@show ineq
#@show all_species
idx_derived = sort([dict[species_derived[i]] for i in eachindex(species_derived)])

idx_integrate_a = filter!(s -> all(s .!= idx_derived), collect(1:N_spec))
const N_integrate = length(idx_integrate_a)
const idx_integrate = SVector{N_integrate}( idx_integrate_a )


#always present
const iH = dict["H"]
const ielec = dict["e-"]
const iH2 = dict["H2"]
const iHp = dict["H+"]

#optional elements
const iHe = any(all_species.=="He") ? dict["He"] : 0
const iC  = any(all_species.=="C" ) ? dict["C" ] : 0
const iO  = any(all_species.=="O" ) ? dict["O" ] : 0
const iSi = any(all_species.=="Si") ? dict["Si"] : 0
const iS  = any(all_species.=="S" ) ? dict["S" ] : 0
const iN  = any(all_species.=="N" ) ? dict["N" ] : 0
const iMg = any(all_species.=="Mg" ) ? dict["Mg" ] : 0
const iFe = any(all_species.=="Fe" ) ? dict["Fe" ] : 0

const iCO  = any(all_species.=="CO" ) ? dict["CO" ] : 0

const iHep = any(all_species.=="He+") ? dict["He+"] : 0
const iCp = any(all_species.=="C+") ? dict["C+"] : 0
const iSip = any(all_species.=="Si+") ? dict["Si+"] : 0

num_a = zeros(Int32, N_reac)
typ_a = zeros(Int32, N_reac)

ir1_a = zeros(Int32, N_reac)
ir2_a = zeros(Int32, N_reac)
ip1_a = zeros(Int32, N_reac)
ip2_a = zeros(Int32, N_reac)
ip3_a = zeros(Int32, N_reac)
ip4_a = zeros(Int32, N_reac)

alpha_a = zeros(N_reac)
beta_a  = zeros(N_reac)
gamma_a = zeros(N_reac)

const iCion = 0
const SupTh = false
#const SupTh = true

for i in 1:N_reac
    num_a[i], t, R, P, alpha_a[i], beta_a[i], gamma_a[i] = read_umist_line(reduced_network[i,:])
    if any(R .== "C+") && any(R .== "e-")
        #alpha_a[i] = 4.67e-12
        #gamma_a[i] = 3.76
        #@show alpha_a[i], beta_a[i], gamma_a[i]
    end
    if any(R .== "H-") && any(R .== "PHOTON")
        #alpha_a[i] = 55.8e-10
        #gamma_a[i] = 3.76
        #@show alpha_a[i], beta_a[i], gamma_a[i]
    end
    if any(R .== "PHOTON") && any(R .== "C") #C photoionization
        iCion = i
        #@show alpha_a[i], gamma_a[i]
    end

    if any(R .== "PHOTON")
        typ_a[i] = 1
    elseif any(R .== "CRP")
        typ_a[i] = 2
    elseif any(R .== "CRPHOT")
        typ_a[i] = 3
    else
        typ_a[i] = 0
    end

    #overwrite if it's a grain-process
    if t == "GR"
        typ_a[i] = 4
    end

    #overwrite if it's an ion-neutral process
    if SupTh == true && t == "IN"
        if any(R .== "C+") && any(R .== "H2") && any(P .== "CH+") && any(P .== "H")
            typ_a[i] = 5
            @show i,R,P
        end
    end

    ir1_a[i], ir2_a[i] = dict[R[1]], dict[R[2]]
    ip1_a[i], ip2_a[i], ip3_a[i], ip4_a[i] = dict[P[1]], dict[P[2]], dict[P[3]], dict[P[4]]
end

const typ = SVector{N_reac}(typ_a); const ir1 = SVector{N_reac}(ir1_a); const ir2 = SVector{N_reac}(ir2_a);
const ip1 = SVector{N_reac}(ip1_a); const ip2 = SVector{N_reac}(ip2_a); const ip3 = SVector{N_reac}(ip3_a);
const ip4 = SVector{N_reac}(ip4_a); const pro4 = any(ip4_a .!= 0); const pro3 = any(ip3_a .!= 0);

const alpha = SVector{N_reac}(alpha_a);
const beta  = SVector{N_reac}(beta_a);
const gamma = SVector{N_reac}(gamma_a);

function setup_conservation_arrays(elem)
    idx_spec = Int32[]
    count_spec = Int32[]
    fac = zeros(Int32, N_spec)

    for i in eachindex(elements)
        for j in eachindex(elements[i])
            #@show i, j, elements[i][j]
            if elements[i][j] == elem
                push!(idx_spec, i)
                push!(count_spec, counts[i][j])
                fac[i] = counts[i][j]
                #@show counts[i][j], all_species[i]
            end
        end
    end
    fac
end


fac_a_H  = setup_conservation_arrays("H"); fac_a_He = setup_conservation_arrays("He");
fac_a_C  = setup_conservation_arrays("C"); fac_a_O  = setup_conservation_arrays("O");
fac_a_Si = setup_conservation_arrays("Si"); fac_a_S  = setup_conservation_arrays("S");
fac_a_N  = setup_conservation_arrays("N"); fac_a_Mg  = setup_conservation_arrays("Mg");
fac_a_Fe  = setup_conservation_arrays("Fe");

const fac_H  = SVector{N_spec}(fac_a_H)
const fac_He = SVector{N_spec}(fac_a_He)
const fac_C  = SVector{N_spec}(fac_a_C)
const fac_O  = SVector{N_spec}(fac_a_O)
const fac_Si = SVector{N_spec}(fac_a_Si)
const fac_S  = SVector{N_spec}(fac_a_S)
const fac_N  = SVector{N_spec}(fac_a_N)
const fac_Mg = SVector{N_spec}(fac_a_Mg)
const fac_Fe = SVector{N_spec}(fac_a_Fe)
const charge = SVector{N_spec}(charge_a);

const XHe = 0.1
const omega = 0.5;

const SEC_PER_YEAR = 31556952.0

#const abC_s = 2.9e-4
#const abO_s = 4.9e-4
#const abN_s = 6.8e-4
#const abSi_s = 3.2e-5

#const abC_s = 1.4e-4
#const abO_s = 3.2e-4
#const abSi_s = 1.5e-5
#const abS_s = 8.0e-6
#const abN_s = 7.5e-5
#const abMg_s = 3.5e-5
#const abFe_s = 2.8e-5
#const abF_s = 2.e-8
#const abNa_s = 2.e-8

const abC_s = 1.0e-4
const abO_s = 3.0e-4
#const abSi_s = 0.0
#const abSi_s = 1.5e-5
const abSi_s = 1.7e-6 #depleted abundance

#const abS_s = 0.0
const abS_s = 1e-6
#const abN_s = 7.6e-5
#const abFe_s = 1.6e-6
#const abMg_s = 1.4e-5
#const abNa_s = 1.3e-6

const abN_s = 0.0
const abFe_s = 0.0
const abMg_s = 0.0
const abNa_s = 0.0
const abF_s = 0.0

const INTEGRATE = true  # 1 = integration, 0 = root-finding

const ξ0 = 1.36e-17

const kdust = 3e-17

const kH2diss = 5.68e-11 * 0.52
const kCOdiss = 2.43e-10 * 0.48

const gamma_H2 = 3.85
const gamma_CO = 3.51

#const grRec = true
const grRec = false
#const ssC = true
const ssC = false

re = Dict(num_a .=> collect(1:N_reac) );
function calc_coeff(coeff, temp, nH, NH, NH2, NCO, NC, xelec, ξp, IUV, Zp)
    Av = NH * 5.35e-22
    for i in 1:N_reac
        if typ[i] == 1
            coeff[i] = IUV * alpha[i] * exp(-gamma[i] * Av)
            #coeff[i] = IUV * exp(-38.0 * Zp) * alpha[i] * exp(-gamma[i] * Av)
        elseif typ[i] == 2
            coeff[i] = alpha[i] * ξp
        elseif typ[i] == 3
            coeff[i] = alpha[i] * ξp * (temp / 300.)^beta[i] * gamma[i] / (1.0 - omega)
        elseif typ[i] == 0
            coeff[i] = alpha[i] * (temp / 300.)^beta[i] * exp(-gamma[i] / temp)
        end
        #overwrite
        if SupTh == true && typ[i] == 5 && NH2 < 4e0
            #=
            m_reduce = 1.714  #m_C+ = 12, m_H2 = 2, 12*2/(12+2) = 1.714
            vA = 3.3 #km/sec
            temp_eff = temp + 40.3825 * m_reduce * vA^2
            coeff[i] = alpha[i] * (temp_eff / 300.)^beta[i] * exp(-gamma[i] / temp_eff)
            if isnan(coeff[i])
                @show "nan!!!", i, coeff[i], temp_eff
            end
            =#
            coeff[i] = alpha[i] * (800. / 300.)^beta[i] * exp(-gamma[i] / 800.)
        end
        if isnan(coeff[i])
            @show i
        end
    end
    #=
    if temp <= 261.
        coeff[re[5643]] = 3.5e-11
    else
        coeff[re[5643]] = 1.77e-11 * exp(178.0/temp)
    end
    lnT = log(temp)
    coeff[re[6162]] = 2.753e-14 * (315614.0/temp)^1.5 * (1.0 + (115188.0/temp)^0.407)^(-2.242)
    coeff[re[6166]] = 1e-11 * temp^(-0.5) * (11.19 - 1.676*lnT - 0.2852*lnT^2 + 0.04433*lnT^3)
    coeff[re[3431]] = 1.4e-9 * (temp/300.)^(-0.5)
    coeff[re[1158]] = 7.0e-8 * (temp/300.)^(-0.5)
    coeff[re[1293]] = 1.08e-7 * (temp/300.)^(-0.5)
    coeff[re[1295]] = 6.02e-8 * (temp/300.)^(-0.5)
    coeff[re[1296]] = 2.58e-7 * (temp/300.)^(-0.5)
    =#
    #photodissociation
    coeff[end]   = IUV * kH2diss * exp(-gamma_H2 * Av) * fH2selfshield(NH2)
    coeff[end-1] = IUV * kCOdiss * exp(-gamma_CO * Av) * fCOselfshield(NCO,NH2)

    #grain processes
    coeff[end-2] = kdust * (temp / 100.)^0.5 * Zp

    psi = 1e20 #large value will turn off grain recombination
    if grRec == true
        if xelec > 1e-20
            psi = 1.7 * (IUV+1e-20) * exp(-1.87*Av) * sqrt(temp) / (nH * xelec)
        end
        psi = psi < 1e20 ? psi : 1e20
    end

    if psi <= 0.0
        println("IUV=", IUV, "  AV=", Av, "  temp=", temp, "  nH=", nH, "  xelec=", xelec)
        error("psi<=0")
    end
    coeff[end-3] = grain_recomb_H( temp, psi) * Zp
    coeff[end-4] = grain_recomb_He(temp, psi) * Zp
    coeff[end-5] = grain_recomb_C( temp, psi) * Zp
    coeff[end-6] = grain_recomb_Si(temp, psi) * Zp

    if ssC == true
        coeff[iCion] *= fCselfshield(NC,NH2)
    end
    #=
    if SupTh == true && NH2 < 4e20
        m_reduce = 1.714  #m_C+ = 12, m_H2 = 2, 12*2/(12+2) = 1.714
        vA = 3.3 #km/sec
        temp_eff = temp + 40.3825 * vA^2
        coeff[iSupTh] = alpha[iSupTh] * (temp_eff / 300.)^beta[iSupTh] * exp(-gamma[iSupTh] / temp_eff)
    end
    =#
    #if fCselfshield(NC,NH2) < 0.1 @show NC, NH2, NCO end
end

function init_abund(abund, Zp)
    calc_abund_derived(abund, Zp)
    abund[iH]=0.4
    abund[iH2]=0.3
end

function calc_abund_derived(abund, Zp)
    if NONEQ
        for i in 1:N_neq
            abund[ineq[i]] = xneq[i]
        end
    end

    sH = sHe = sC = sO = sSi = sS = sN = sMg = sFe = selec = 0.0
    for i in 1:N_spec
        (iH  > 0) && (i != iH)  ? (sH  += abund[i] * fac_H[i])  : nothing
        (iHe > 0) && (i != iHe) ? (sHe += abund[i] * fac_He[i]) : nothing
        (iC  > 0) && (i != iC)  ? (sC  += abund[i] * fac_C[i])  : nothing
        (iO  > 0) && (i != iO)  ? (sO  += abund[i] * fac_O[i])  : nothing
        (iSi > 0) && (i != iSi) ? (sSi += abund[i] * fac_Si[i]) : nothing
        (iS  > 0) && (i != iS)  ? (sS  += abund[i] * fac_S[i])  : nothing
        (iN  > 0) && (i != iN)  ? (sN  += abund[i] * fac_N[i])  : nothing
        (iMg > 0) && (i != iMg) ? (sMg += abund[i] * fac_Mg[i]) : nothing
        (iFe > 0) && (i != iFe) ? (sFe += abund[i] * fac_Fe[i]) : nothing
        (ielec  > 0) && (i != ielec) ? (selec += abund[i] * charge[i]) : nothing
    end

    iH  > 0 ? abund[iH]  = 1.0 - sH : nothing
    iHe > 0 ? abund[iHe] = XHe - sHe : nothing
    iC  > 0 ? abund[iC]  = abC_s * Zp - sC : nothing
    iO  > 0 ? abund[iO]  = abO_s * Zp - sO : nothing
    iS  > 0 ? abund[iS]  = abS_s * Zp - sS : nothing
    iSi > 0 ? abund[iSi] = abSi_s * Zp - sSi : nothing
    iN  > 0 ? abund[iN]  = abN_s * Zp - sN : nothing
    iMg > 0 ? abund[iMg] = abMg_s * Zp - sMg : nothing
    iFe > 0 ? abund[iFe] = abFe_s * Zp - sFe : nothing
    ielec > 0 ? abund[ielec] = selec : nothing

end


function solve_equilibrium_abundances(abund, dtime, nH, temp, NH, NH2, NCO, NC, ξ, IUV, Zp)
#function solve_equilibrium_abundances(abund, nH, temp, NH, NH2, NCO, NC, ξ, IUV, Zp)

    ξp = ξ / ξ0;

    params = [nH, temp, NH, NH2, NCO, NC, ξ, IUV, Zp]
    #params = SVector{6,Float64}(nH, temp, NH, IUV, ξp, Zp)


    #initial conditions
    #abund = zeros(N_spec)
    #init_abund(abund, Zp)
    #@show abund

    abund_dot = zeros(N_spec)
    abund_final = zeros(N_spec)
    abund_con = zeros(N_integrate)

    #coeff = zeros(N_reac+7) #+1 for H2 formation on dust
    coeff = zeros(N_reac) #+1 for H2 formation on dust
    reac_rates = zeros(N_reac)

    if INTEGRATE == true
        function calc_abund_dot_closure(abund_dot_int, abund_int, ppp, t)

            abund_dot .= 0.0 #in-place for better performance

            abund[idx_integrate] .= abund_int
            calc_abund_derived(abund, Zp)

            calc_coeff(coeff, temp, nH, NH, NH2, NCO, NC, abund[ielec], ξp, IUV, Zp) #zero allocation

            for i in 1:N_reac
                if typ[i] == 0
                    reac_rates[i] = coeff[i] * abund[ir1[i]] * abund[ir2[i]] * nH  #2-body process
                elseif typ[i] == 3
                    reac_rates[i] = coeff[i] * abund[ir1[i]] * abund[iH2]  #CR-induced photon process
                elseif typ[i] == 4
                    reac_rates[i] = coeff[i] * abund[ir1[i]] * nH  #grain processes
                else
                    #typ[i] == 1 or 2
                    reac_rates[i] = coeff[i] * abund[ir1[i]]  #CR or photo process
                end
            end #zero allocation

            for i in 1:N_reac

                #destruction
                if ir1[i] != 0
                    abund_dot[ir1[i]] -= reac_rates[i]
                end
                if ir2[i] != 0
                    abund_dot[ir2[i]] -= reac_rates[i]
                end

                #formation
                if ip1[i] != 0
                    abund_dot[ip1[i]] += reac_rates[i]
                end
                if ip2[i] != 0
                    abund_dot[ip2[i]] += reac_rates[i]
                end
                if pro3 == true
                    if ip3[i] != 0
                        abund_dot[ip3[i]] += reac_rates[i]
                    end
                end
                if pro4 == true
                    if ip4[i] != 0
                        abund_dot[ip4[i]] += reac_rates[i]
                    end
                end
            end #zero allocation
            abund_dot_int .= abund_dot[idx_integrate]
        end
        #tend = SEC_PER_YEAR * 1e10 / Zp / nH
        tend = SEC_PER_YEAR * dtime
        tspan = (0.0, tend)
        abund_con .= abund[idx_integrate] #can't use idx_integrate
        #prob = ODEProblem(calc_abund_dot_wrap, abund_con, tspan, params)
        prob = ODEProblem(calc_abund_dot_closure, abund_con, tspan)
        #sol = solve(prob, Rosenbrock23(autodiff=false), reltol=1e-9, abstol=1e-9);
        #sol = solve(prob, CVODE_BDF(), reltol=1e-7, abstol=1e-7);
        sol = solve(prob, reltol=1e-7, abstol=1e-7);
        #sol = solve(prob, reltol=1e-13, abstol=1e-13);
        #sol = solve(prob, isoutofdomain=(u,p,t)->any(x->x<0,u), reltol=1e-9, abstol=1e-9);
        #cb = ManifoldProjection(g)
        #sol = solve(prob, Rosenbrock23(autodiff=false), save_everystep=false, callback=cb)
        abund_final[idx_integrate] .= sol.u[end]
        calc_abund_derived(abund_final, Zp)
    else
        #===
        function f!(xdot, x)
            calc_abund_dot(xdot, x, params, 0)
            xdot .*= 1e15
            xdot[ielec] = sum( x .* charge )
            xdot[iH]  = 1.0 - sum( x[idx_spec_H]  .* count_spec_H  )
            #@show x
            #@show xdot
        end
        sol = nlsolve(f!, abund, autodiff = :forward, ftol=1e-10, iterations=1000)
        @show sol.x_converged
        abund_final .= sol.zero
        ===#
    end
    abund .= abund_final
    return abund_final, reac_rates
end

end
