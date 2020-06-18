species_noneq = []
#species_noneq = ["H2"]
#species_noneq = ["H2", "H+"]
const N_neq = length(species_noneq)
const NONEQ = N_neq > 0

const file_path = "/ptmp/huchiayu/snapshots/tallbox/SFSNPI_N1e6_gS10H250dS40H250_soft2_SFMJ1_eff0p5_60Myr_Lsh100"
const Zp = 1.0
const ShieldingLength = 0.1
include("treecol.jl")

#i = 640
#abund_all, ga, X, NH, NH2, NCO, tree = solve_chem_all_particles(i);



#const Zp = parse(T, ARGS[2])

@show file_path
@show Zp

#snaps = collect(800:-10:150)
snaps = collect(150:10:1000)
for i in snaps
    abund_all, ga, X, NH, NH2, NCO, tree = solve_chem_all_particles(i, file_path, Zp);
end
0



