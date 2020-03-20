using Interpolations

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
    return exp( -1.6e-17 * NC ) * exp(-rH2) / (1.0 + rH2)
end
