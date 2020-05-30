module OctTree
using StaticArrays
using Healpix
const NSIDE = 1
const NPIX = 12*NSIDE^2


export Node, TreeGather
export buildtree, get_scatter_ngb_tree, get_gather_ngb_tree, treewalk, nearest
#export plot_quadtree, plot_treewalk, mycircle, plot_circles_scatter_ngbs
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
	pos_c::SVector{N,T} #center of mass
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
			node.n.pos_c = p.pos
            node.p = p
			#@show "empty leaf!!!", node.n.mass, p.mass, node.p.pos, node.n.pos_c
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
			masssum = node.n.mass + p.mass
			inv_masssum = 1.0 / masssum
			node.n.pos_c = (node.n.mass .* node.n.pos_c .+ p.mass .* p.pos) .* inv_masssum
			#@show "non-empty leaf!!!", node.n.mass, p.mass, sum((node.n.pos_c - p.pos).^2)^0.5 / node.length[1]
            node.n.mass = masssum #this line has to be put after insertpart!, otherwise we're inserting back the wrong node.mass!!!
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
		masssum = node.n.mass + p.mass
		inv_masssum = 1.0 / masssum
		node.n.pos_c = (node.n.mass .* node.n.pos_c .+ p.mass .* p.pos) .* inv_masssum
		#@show "nonleaf!!!", node.n.mass, p.mass, sum((node.n.pos_c - p.pos).^2)^0.5 / node.length[1]
        node.n.mass = masssum
		node.n.mass_H2 += p.mass_H2
		node.n.mass_CO += p.mass_CO
        insertpart!(p, node.child[getChildIndex(p.pos,node)])
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

#when changing the argument of treewalk, remember to also change the recursive treewalk call inside the function body!
function treewalk(ga::TreeGather{T}, p::SVector{N,T}, node::Node{N,T}, openingangle::T, shieldinglength::T, boxsizes::SVector{N,T}) where {N,T}
    if isLeaf(node)
        #println("in a leaf node")
        if node.p == nothing
            #println("empty leaf")
        else
            #println("nonempty leaf")
            ga.mass += node.n.mass
            dx = nearest.(node.p.pos - p, boxsizes)
			dist2 = sum(dx.^2)
			#@show dist2, node.p.pos, p
			if 0.0 < dist2 < shieldinglength^2  #exclude self contribution (dist2=0)
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
            dist2 = get_distance2(node.child[i].n.pos_c, p, boxsizes, true)
            if dist2 > (node.child[i].length[1] / openingangle)^2
                #println("skip node ", i)
				if dist2 < shieldinglength^2
	                ga.mass += node.child[i].n.mass
	                dx = nearest.(node.child[i].n.pos_c - p, boxsizes)
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
            treewalk(ga, p, node.child[i], openingangle, shieldinglength, boxsizes)
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

end #module OctTree
