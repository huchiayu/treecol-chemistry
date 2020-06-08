push!(LOAD_PATH, pwd())
#import Pkg; Pkg.add("HDF5"); Pkg.add("StaticArrays"); Pkg.add("PyPlot");
#Pkg.add("DifferentialEquations"); Pkg.add("Parameters"); Pkg.add("Sundials");
using OctTree

using HDF5
using StaticArrays
using Statistics
using LinearAlgebra
using .Threads
#using Serialization
using Random
using BenchmarkTools

#import Plots #if needed, it has to be imported before PyPlot otherwise it'll crash
using PyPlot


const BOXSIZE_X = 1.0
const BOXSIZE_Y = 1.0
const BOXSIZE_Z = 1.0

const ANGLE = 0.7
const ShieldingLength = 0.1

const N=3
const T=Float64


function plot_quadtree(node::Node{N,T}, ix, iy) where {N,T}
    if N<2 println("N must be >= 2")
        return
    end
    xmin = node.center[ix] - 0.5*node.length[ix]
    xmax = node.center[ix] + 0.5*node.length[ix]
    ymin = node.center[iy] - 0.5*node.length[iy]
    ymax = node.center[iy] + 0.5*node.length[iy]
    color="grey"
    plot([xmin,xmin],[ymin,ymax], c=color)
    plot([xmin,xmax],[ymin,ymin], c=color)
    plot([xmax,xmax],[ymin,ymax], c=color)
    plot([xmin,xmax],[ymax,ymax], c=color)
    if node.child != nothing
        for i in 1:2^N
            plot_quadtree(node.child[i], ix, iy)
        end
    end
end

function plot_circles_scatter_ngbs(X::Vector{SVector{N,T}}, hsml::Vector{T}, boxsizes::SVector{N,T}) where {N,T}
    for i in eachindex(X)
        mycircle(hsml[i], X[i][1], X[i][2], 0.5*boxsizes[1], 0.5*boxsizes[2], boxsizes[1], boxsizes[2], "green")
    end
end

function mycircle(r, xc, yc, x0, y0, boxsizeX, boxsizeY, color)
    #xc=1.5;yc=2.5;r=0.5
    x = collect(-r : 0.0001*r : r)
    y = sqrt.( r.^2 .- x.^2 )
    xx = x.+xc
    yy = y.+yc
    #x0 = 0.5*boxHalf_X # coordinates for the center of box
    #y0 = 0.5*boxHalf_Y # coordinates for the center of box
    xx[xx .> x0+boxsizeX*0.5] .-= boxsizeX
    xx[xx .< x0-boxsizeX*0.5] .+= boxsizeX
    yy[yy .> y0+boxsizeY*0.5] .-= boxsizeY
    yy[yy .< y0-boxsizeY*0.5] .+= boxsizeY
    scatter(xx,yy,marker=".",c=color,s=0.03)
    xx = x.+xc
    yy = -y.+yc
    xx[xx .> x0+boxsizeX*0.5] .-= boxsizeX
    xx[xx .< x0-boxsizeX*0.5] .+= boxsizeX
    yy[yy .> y0+boxsizeY*0.5] .-= boxsizeY
    yy[yy .< y0-boxsizeY*0.5] .+= boxsizeY
    scatter(xx,yy,marker=".",c=color,s=0.03)
end

function plot_treewalk(ga::TreeGather{T},ix,iy) where {T}
    for i in eachindex(ga.nodecenters)
        if ga.nodelengths[i][1] == 0
            #plot(ga.nodecenters[i][ix], ga.nodecenters[i][iy], ".", c="red")
        else
            xmin = ga.nodecenters[i][ix] - 0.5*ga.nodelengths[i][ix]
            xmax = ga.nodecenters[i][ix] + 0.5*ga.nodelengths[i][ix]
            ymin = ga.nodecenters[i][iy] - 0.5*ga.nodelengths[i][iy]
            ymax = ga.nodecenters[i][iy] + 0.5*ga.nodelengths[i][iy]
            color="green"
            plot([xmin,xmin],[ymin,ymax], c=color)
            plot([xmin,xmax],[ymin,ymin], c=color)
            plot([xmax,xmax],[ymin,ymax], c=color)
            plot([xmin,xmax],[ymax,ymax], c=color)
        end
    end
end

if N==3
    const boxsizes = SVector{N,T}(BOXSIZE_X, BOXSIZE_Y, BOXSIZE_Z)
elseif N==2
    const boxsizes = SVector{N,T}(BOXSIZE_X, BOXSIZE_Y)
elseif N==1
    const boxsizes = SVector{N,T}(BOXSIZE_X)
end

function treecol_test(X::Vector{SVector{N,T}}) where {N,T}
    #boxsizes = SVector{N,T}(1.0,1.0,1.0) .* 2
    hsml0 = 0.1*BOXSIZE_X
    Npart = length(X)
    hsml = ones(T,Npart) .* hsml0
    hsml .*= (1.5 .- rand(Npart))
    mass = ones(Npart)
    mass_H2 = ones(Npart)
    mass_CO = ones(Npart)
    @time tree = buildtree(X, hsml, mass, mass_H2, mass_CO, boxsizes);


    ga_out = Vector{TreeGather{Float64}}(undef,Npart)

    @time for i in 1:Npart
        ga = TreeGather{T}()
        treewalk(ga, X[i], tree, ANGLE, ShieldingLength, boxsizes)
        ga_out[i] = ga
    end


    return tree, ga_out
end

#=
using Healpix
X = [@SVector rand(N) for _ in 1:100000]
push!(X,SVector(0.5,0.5,0.5).*BOXSIZE_X) #add a particle at the box center
tree, ga_out = treecol_test(X);

const fac_col = 8.8674721e23
m = Map{Float64, RingOrder}(NSIDE);
m.pixels[:] = log10.(ga_out[end].column_all .* fac_col) ;
#Plots.plot(m, clim=(20.5,24))
Plots.plot(m)
=#

function ngb_test(p, hsml0)
    #boxsizes = SVector{N,T}(1.0,1.0,1.0)
    #idx_ngbs = get_scatter_ngb_tree(p, tree, boxsizes)
    idx_ngbs = get_gather_ngb_tree(p, hsml0, tree, boxsizes)
    @show length(idx_ngbs)
    ix,iy=1,2
    xngb = getindex.(X[idx_ngbs], ix)
    yngb = getindex.(X[idx_ngbs], iy)
    fig = figure("tree",figsize=(8,8))
    #plot_quadtree(tree, ix, iy)  #costly!!!
    scatter(getindex.(X,ix), getindex.(X,iy), c="blue", marker="o", s=3)
    scatter( xngb, yngb , marker="o", color="red", s=10)
    #@show xngb, yngb
    mycircle(hsml0, p[ix], p[iy], 0.5*boxsizes[ix], 0.5*boxsizes[iy], boxsizes[ix], boxsizes[iy], "red")
    #plot_circles_scatter_ngbs(X[idx_ngbs], hsml[idx_ngbs], boxsizes)
    title("neighbor finder with quadtree")
    xlabel("x")
    ylabel("y")
end
#p=SVector{N,T}(0.,0.5,0.5)
#ngb_test(p,hsml0)

function loop_all_particles_ngbs(X::Vector{SVector{N,T}}, hsml0::T) where {N,T}

    Npart = length(X)
    hsml = ones(T,Npart) .* hsml0
    #hsml .*= (1.5 .- rand(Npart))
    mass = ones(Npart)
    mass_H2 = ones(Npart)
    mass_CO = ones(Npart)
    #tree = buildtree(X, hsml, mass, mass_H2, mass_CO, boxsizes);
    topnode_length = SVector{3}(4.,4.,4.)
    center = topnode_length .* 0.5
    tree = buildtree(X, hsml, mass, mass_H2, mass_CO, center, topnode_length);

    Nngbs = zeros(Int64,length(X))
    for i in eachindex(X)
    #for i in 1:1
        idx_ngbs = get_gather_ngb_tree(X[i], hsml0, tree, boxsizes)
        #idx_ngbs = get_scatter_ngb_tree(X[i], tree, boxsizes)
        Nngbs[i] = length(idx_ngbs)
    end
    Nngbs,tree
end

X = [@SVector rand(N) for _ in 1:100000]
#hsml0 = 0.04
fac_geo = 4*pi/3
Nngb0 = 32
hsml0 = BOXSIZE_X * (Nngb0/(fac_geo*length(X)))^(1/N)
@time Nngbs,tree = loop_all_particles_ngbs(X,hsml0)
#@btime Nngbs = loop_all_particles_ngbs($X,$hsml0)


using NearestNeighbors
function loop_all_particles_ngbs_kdtree(X, hsml0)
    data = zeros(N,length(X))
    for i in eachindex(X) data[:,i] = X[i] end
    kdtree = BallTree(data,PeriodicEuclidean(SVector(1.0, 1.0, 1.0)))
    Nngbs = zeros(Int64,length(X))
    for i in eachindex(X)
        idx_ngbs = inrange(kdtree, data[:,i], hsml0, false)
        Nngbs[i] = length(idx_ngbs)
    end
    Nngbs
end

@time NngbsKD = loop_all_particles_ngbs_kdtree(X,hsml0)
#@btime NngbsKD = loop_all_particles_ngbs_kdtree($X,$hsml0)

@show sum(Nngbs .!= NngbsKD)

#println("c.o.m. of the top node = ", sum(X)/length(X) .== tree.n.pos_c)
#@show sum(X)/length(X) .≈ tree.n.pos_c
