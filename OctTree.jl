module OctTree
using StaticArrays
using PyPlot
using Healpix
const NSIDE = 1
const NPIX = 12*NSIDE^2


export Node, TreeGather
export buildtree, get_scatter_ngb_tree, get_gather_ngb_tree, treewalk, nearest
export plot_quadtree, plot_treewalk
export NSIDE, NPIX

mutable struct PartData{N,T}
	pos::SVector{N,T}
	idx::Int64
	hsml::T
	mass::T
	mass_H2::T
	mass_CO::T
end

PartData{N,T}() where {N,T} = PartData{N,T}(zero(SVector{N,T}),0,0,0,0,0)

mutable struct NodeData{N,T} #auxiliary data carried by treenodes
	com::SVector{N,T} #center of mass
	max_hsml::T  #for scatter ngb search
	mass::T
	mass_H2::T
	mass_CO::T
end

NodeData{N,T}() where {N,T} = NodeData{N,T}(zero(SVector{N,T}),0,0,0,0)

mutable struct Node{N,T}
    center::SVector{N,T}
    length::SVector{N,T}
	#com::SVector{N,T} #center of mass
    #part::Union{SVector{N,T}, Nothing}
    #part_idx::Int64
    #max_hsml::T #for scatter ngb search
    #mass::T
	n::NodeData{N,T}
	p::Union{PartData{N,T}, Nothing}
    child::Union{Array{Node{N,T},1}, Nothing}
end

Node{N,T}(center::SVector{N,T}, length::SVector{N,T}) where {N,T} =
Node{N,T}(center::SVector{N,T}, length::SVector{N,T}, NodeData{N,T}(), nothing, nothing)

function isLeaf(node::Node{N,T}) where{N,T}
    return node.child == nothing ? true : false
end

function getChildIndex(pos::SVector{N,T}, node::Node{N,T}) where {N,T}
    idx = 0
    for i in 1:N
        if pos[i] > node.center[i]
            idx |= 2^(i-1)
        end
        #println(i, 2^(i-1), idx)
    end
    return idx+1
end

function getoffset(i,N)
#for i in 1:2^N
    a = bitstring(i-1)
    #offset = SVector(2*parse(Int, a[end])-1, 2*parse(Int, a[end-1])-1, 2*parse(Int, a[end-2])-1)
    offset = Int[]
    for j in 1:N
        push!(offset,2*parse(Int, a[end-(j-1)])-1)
    end
    #println(offset)
    return offset
end

#function insertpart!(p::SVector{N,T}, p_idx::Int64, node::Node{N,T}) where {N,T}
#function insertpart!(p::PartData{N,T}, p_idx::Int64, hsml::T, mass::T, node::Node{N,T}) where {N,T}
function insertpart!(p::PartData{N,T}, node::Node{N,T}) where {N,T}
    if isLeaf(node)
        if node.p == nothing
            #println("at an empty leaf, insert particle and we're done")
            node.p = p
			node.n.max_hsml = p.hsml
			node.n.mass = p.mass
			node.n.mass_H2 = p.mass_H2
			node.n.mass_CO = p.mass_CO
	    	#node.part_idx = p_idx
            #node.max_hsml = p.hsml
            #node.mass = mass
        else
            #println("at a leaf with a preexisting particle, so we split it")
            node.child = Array{Node{N,T},1}(undef,2^N)
            for i in 1:2^N
                #println(getoffset(i,N))
                childcenter = node.center + getoffset(i,N) * 0.25 .* node.length
                #println("childcenter=", childcenter, "  getoffset(i,N)=", getoffset(i,N))
                #node.child[i] = Node{N,T}(childcenter, 0.5*node.length, nothing, 0, 0, 0, nothing)
				node.child[i] = Node{N,T}(childcenter, 0.5*node.length)
                #println("center = ", childcenter, "  edge_min = ", childcenter - 0.25*node.length,
                #        "  edge_max = ", childcenter + 0.25*node.length)
            end
            if p.hsml > node.n.max_hsml
                node.n.max_hsml = p.hsml
            end
            #println("insert back the preexsisting particle...")
            #insertpart!(node.part, node.part_idx, node.max_hsml, node.mass, node.child[getChildIndex(node.part, node)])
			insertpart!(node.p, node.child[getChildIndex(node.p.pos, node)])
            node.n.mass += p.mass #this line has to be put after the above line, otherwise we're inserting back the wrong node.mass!!!
			node.n.mass_H2 += p.mass_H2
			node.n.mass_CO += p.mass_CO
            #println("insert the new particle...")
            insertpart!(p, node.child[getChildIndex(p.pos, node)])
        end
    else
        #println("open node")
        if p.hsml > node.n.max_hsml
            node.n.max_hsml = p.hsml
        end
        node.n.mass += p.mass
		node.n.mass_H2 += p.mass_H2
		node.n.mass_CO += p.mass_CO
        insertpart!(p, node.child[getChildIndex(p.pos,node)])
    end
end

function plot_quadtree(node::Node{N,T}, ix, iy) where {N,T}
    #println("center=", node.center)
    #println("length=", node.length)
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

#const boxHalf_X = 0.5

#function nearest(x::T) where {T}
function nearest(x::T, boxsize::T) where {T}
    return (x > 0.5*boxsize) ? (x - boxsize) : ( (x < -0.5*boxsize) ? (x + boxsize) : x )
end

function nearest_plain(x::T) where {T}
    if x < -boxHalf_X
        x = x + boxHalf_X*2
    else
        x = x
    end
    if x > boxHalf_X
        x = x - boxHalf_X*2
    else
        x = x
    end
    return x
end

#function get_distance2(a::SVector{N,T}, b::SVector{N,T}, periodic::Bool) where {N,T}
function get_distance2(a::SVector{N,T}, b::SVector{N,T}, boxsizes::SVector{N,T}, periodic::Bool) where {N,T}
    c = a - b
    if periodic
        c = nearest.(c, boxsizes)
    end
    return sum(c .^ 2)
    #===
    for i in 1:N
        c = a[i] - b[i]
        if periodic
            c = nearest(c, 0.5*boxsize[i])
        end
        sum += c^2
    end
    ===#
end

res = Healpix.Resolution(NSIDE)

function vec2pix(p)
    theta,phi = vec2ang(p[1],p[2],p[3])
    ang2pixRing(res, theta,phi)
end
#const ANGLE = 0.7
const ShieldingLength = 0.1

mutable struct TreeGather{T}
    mass::T
    column_all::Vector{T}
    column_H2::Vector{T}
    column_CO::Vector{T}

    #nodecenters::Vector{Vector{T}} #for debugging/visualization
    #nodelengths::Vector{Vector{T}} #for debugging/visualization
end

#TreeGather{T}() where {T} = TreeGather{T}(0,zeros(NSIDE^2*12),zeros(NSIDE^2*12),zeros(NSIDE^2*12),[],[])
TreeGather{T}() where {T} = TreeGather{T}(0,zeros(NPIX),zeros(NPIX),zeros(NPIX))

#function treewalk(ga::TreeGather{T}, p::SVector{N,T}, node::Node{N,T}, boxsizes::SVector{N,T}) where {N,T}
function treewalk(ga::TreeGather{T}, p::SVector{N,T}, node::Node{N,T}, openingangle::T, boxsizes::SVector{N,T}) where {N,T}
    if isLeaf(node)
        #println("in a leaf node")
        if node.p == nothing
            #println("empty leaf")
        else
            #println("nonempty leaf")
            ga.mass += node.n.mass
            dx = nearest.(node.center - p, boxsizes)
			dist2 = sum(dx.^2)
			if dist2 < ShieldingLength^2
	            ipix = vec2pix(dx)
	            #ga.column_all[ipix] += node.mass
	            area = 4*pi/NPIX * dist2
				ga.column_all[ipix] += node.n.mass / area
				ga.column_H2[ipix] += node.n.mass_H2 / area
	            ga.column_CO[ipix] += node.n.mass_CO / area
	            #push!(ga.nodecenters, node.p.pos)
	            #push!(ga.nodelengths, zeros(N))
			end
            #@show "it's a particle..." node.center, node.length
            #dist2 = get_distance2(node.part, p, boxsizes, true)
        end
    else
        #println("This is a node... check its children")
        for i in 1:2^N
            dist2 = get_distance2(node.child[i].center, p, boxsizes, true)
            if dist2 > (node.child[i].length[1] / openingangle)^2
                #println("skip node ", i)
				if dist2 < ShieldingLength^2
	                ga.mass += node.child[i].n.mass
	                dx = nearest.(node.child[i].center - p, boxsizes)
	                ipix = vec2pix(dx)
	                #ga.column_all[ipix] += node.child[i].mass
	                area = 4*pi/NPIX * dist2
	                ga.column_all[ipix] += node.child[i].n.mass / area
					ga.column_H2[ipix] += node.child[i].n.mass_H2 / area
					ga.column_CO[ipix] += node.child[i].n.mass_CO / area
	                #push!(ga.nodecenters, node.child[i].center)
	                #push!(ga.nodelengths, node.child[i].length)
				end
                #@show "use this node!" node.child[i].center, node.child[i].length
                @goto escape_label
            end
            #println("open this node")
            #@show "!!!", m
            treewalk(ga, p, node.child[i], openingangle, boxsizes)
            @label escape_label
        end
    end
end


function get_ptl_in_box_gather(radius::T, p::SVector{N,T}, node::Node{N,T}, boxsizes::SVector{N,T}, idx_ngbs::Vector{Int64}) where {N,T}
    if isLeaf(node)
        #println("in a leaf node")
        if node.p == nothing
            #println("empty leaf")
        else
            #println("nonempty leaf")
	    dist2 = get_distance2(node.p.pos, p, boxsizes, true)
            if dist2 < radius^2
                #println("push ", node.part)
		#push!(x_ngbs, node.p.pos)
		#push!(dist, sqrt(dist2))
		#push!(idx_ngbs, node.p.idx)
		push!(idx_ngbs, node.p.idx)
            end
        end
    else
        #println("This is a node... check its children")
        for i in 1:2^N
            dist2 = get_distance2(node.child[i].center, p, boxsizes, true)
            for j in 1:N
                if abs(nearest(node.child[i].center[j] - p[j], boxsizes[j])) > 0.5*node.child[i].length[j] + radius
                    #println("skip node ", i)
                    @goto escape_label
                end
            end
            #println("open this node")
            get_ptl_in_box_gather(radius, p, node.child[i], boxsizes, idx_ngbs)
            @label escape_label
        end
    end
end

function get_ptl_in_box_scatter(p::SVector{N,T}, node::Node{N,T}, boxsizes::SVector{N,T}, idx_ngbs::Vector{Int64}) where {N,T}
    if isLeaf(node)
        #println("in a leaf node")
        if node.p == nothing
            #println("empty leaf")
        else
            #println("nonempty leaf")
	    dist2 = get_distance2(node.p.pos, p, boxsizes, true)
            if dist2 < node.n.max_hsml^2
                #println("push ", node.part)
		#push!(x_ngbs, node.part)
		#push!(dist, sqrt(dist2))
		push!(idx_ngbs, node.p.idx)
            end
        end
    else
        #println("This is a node... check its children")
        for i in 1:2^N
            dist2 = get_distance2(node.child[i].center, p, boxsizes, true)
            for j in 1:N
                if abs(nearest(node.child[i].center[j] - p[j], boxsizes[j])) > 0.5*node.child[i].length[j] + node.child[i].n.max_hsml
                    #println("skip node ", i)
                    @goto escape_label
                end
            end
            #println("open this node")
            get_ptl_in_box_scatter(p, node.child[i], boxsizes, idx_ngbs)
            @label escape_label
        end
    end
end

function get_gather_ngb_tree(x::SVector{N,T}, h::T, node::Node{N,T}, boxsizes::SVector{N,T}) where {N,T}
    idx_ngbs = Int64[]
    get_ptl_in_box_gather(h, x, node, boxsizes, idx_ngbs)
    return idx_ngbs
end

function get_scatter_ngb_tree(x::SVector{N,T}, node::Node{N,T}, boxsizes::SVector{N,T}) where {N,T}
    idx_ngbs = Int64[]
    get_ptl_in_box_scatter(x, node, boxsizes, idx_ngbs)
    return idx_ngbs
end

function buildtree(X::Vector{SVector{N,T}}, hsml::Vector{T},
	mass::Vector{T}, mass_H2::Vector{T}, mass_CO::Vector{T}, boxsizes::SVector{N,T}) where {N,T}
    #N_gas::Int32 = sizeof(X)

    #construct the tree for ngb search
    #tree = Node{N,T}(0.5*boxsizes, boxsizes, nothing, 0, 0, 0, nothing)
	tree = Node{N,T}(0.5*boxsizes, boxsizes)

    for i in eachindex(X)
		#insertpart!(SVector(X[i]), i, hsml[i], mass[i], tree)
		part = PartData(SVector(X[i]), i, hsml[i], mass[i], mass_H2[i], mass_CO[i])
		insertpart!(part, tree)
    end
    return tree
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

end #module OctTree

#===
using PyPlot
const BOXSIZE_X = BOXSIZE_Y = BOXSIZE_Z = 1.0

function test(hsml0::T, p::SVector{N,T}, Npart::Int64) where {N,T}

if N==3
    boxsizes = SVector{N,T}(BOXSIZE_X, BOXSIZE_Y, BOXSIZE_Z)
elseif N==2
    boxsizes = SVector{N,T}(BOXSIZE_X, BOXSIZE_Y)
elseif N==1
    boxsizes = SVector{N,T}(BOXSIZE_X)
end

#N=3; T=Float64;
X = [@SVector rand(N) for _ in 1:Npart]

hsml = ones(T,Npart) .* hsml0
hsml .*= (1.5 .- rand(Npart))
tree = buildtree(X, hsml, boxsizes);

#x_ngbs = SVector{N,T}[]
idx_ngbs = get_scatter_ngb_tree(p, tree, boxsizes)
#idx_ngbs = get_gather_ngb_tree(p, hsml0, tree, boxsizes)
X[idx_ngbs];
ix=1
iy=2
xngb = getindex.(X[idx_ngbs], ix)
yngb = getindex.(X[idx_ngbs], iy)

fig = figure("tree",figsize=(8,8))
plot_quadtree(tree, ix, iy)
scatter(getindex.(X,ix), getindex.(X,iy), c="blue", marker="o", s=3)
scatter( xngb, yngb , marker="o", color="red", s=10)
#@show xngb, yngb

mycircle(hsml0, p[ix], p[iy], 0.5*boxsizes[ix], 0.5*boxsizes[iy], boxsizes[ix], boxsizes[iy], "red")
plot_circles_scatter_ngbs(X[idx_ngbs], hsml[idx_ngbs], boxsizes)
#title("neighbor finder with quadtree")
#xlabel("x")
#ylabel("y")
end

p = SVector(0.0, 0.0, 0.0) .+ 0.5
#p = SVector(0.5, 0.5)
radius = 0.2
test(radius, p, 30)
===#


#===
using PyPlot

NN=2; TT=Float64;
Npart = 200
x = (@SMatrix rand(NN,Npart)) - 0.5
#lx = 1.0; ly = 1.0; lz = 1.0;
boxsize = 1.0
node = Node{NN,TT}(zeros(SVector{NN})[1:NN], zeros(SVector{NN})[1:NN] .+ boxsize, nothing, 0, nothing)
for i in 1:Npart
    insertpart!(x[:,i], i, node)
end

x_ngbs = SVector{NN,TT}[]
dist = Float64[]
idx_ngbs = Int64[]
#p = SVector(-0.45, 0.45, 0.45)
#p = SVector(0.0, 0.45)
p = x[:,1]
radius = 0.2
get_ptl_in_box(radius, p, node, x_ngbs, dist, idx_ngbs)

ix=1
iy=2
xngb = x_ngbs[ix:NN:end]
yngb = x_ngbs[iy:NN:end]

fig = figure("tree",figsize=(10,9))
plot_quadtree(node, ix, iy)
scatter(x[ix,:], x[iy,:], c="blue", s=20.)
#scatter(xngb, yngb, c="red", s=20.)
scatter( getindex.(x_ngbs, 1), getindex.(x_ngbs, 2) , marker=".", color="red")
mycircle(radius, p[ix], p[iy], "red")
title("neighbor finder with quadtree")
xlabel("x")
ylabel("y")
===#
