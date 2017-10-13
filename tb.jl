using PyPlot

#= --------------------------------------------------------- =#
mutable struct TBmodel{T<:Real}
    pos::Array{T,2} 			 # position of elements in unit cell
	periodicity::Vector{Bool}	 # periodicity along latx,laty,latz (boolean = yes/no)
    neighbors::Vector{Vector{S}} where S<:Integer # neighboring indices
	hopfun::Function 			 # hopping function (takes two positions)
#    lat::Vector{Vector{T}}       # lattice vectors
end
TBmodel() = TBmodel( Array{Float64}(0,2),
					 Vector{Bool}(0),
					 Vector{Vector{Int}}(0), 
					 (x,y) -> zero(eltype(x)) )


#= --------------------------------------------------------- =#
function interDistance{T<:Real}(pos_i::Vector{T},pos_j::Vector{T},periodicity::Vector{Bool})

	# not periodic in any dimension; regular flat metric
	if iszero(periodicity) 
		return sqrt( sum( (pos_i - pos_j).^2 ) )

	# periodic in some directions; use distance on a flat n-torus,
	# with n = sum(periodicity)
	else 
		s = zero(eltype(pos_i))
		rij = abs.(pos_i - pos_j)
		for xyz in eachindex(pos_i)
			if periodicity[xyz] # dimension is periodic
				s += min(rij[xyz],1-rij[xyz])^2
			else			    # dimension is not periodic
				s += rij[xyz]^2
			end
		end
		return sqrt(s)
	end

end

interDistance{T<:Real}(pos_i::Vector{T},pos_j::Vector{T}) = 
	interDistance(pos_i,pos_j,[false,false,false])


#= --------------------------------------------------------- =#
function setHopping!(tb::TBmodel)
	n = size(tb.pos,1)
	H = spzeros(n,n)

	#= fill upper/lower triangle with hopping amplitude by eval-
	uating the function hopfun as a function of inter-positions  =#
	for i = 1:n
		for j = (i+1):n
			t = tb.hopfun(tb.pos[i,:], tb.pos[j,:])
			if t != 0.0
				H[i,j] = t  # upper triangular part
				H[j,i] = t' # lower triangular part
			end
		end
	end
	# ---- diagonal of H should be computed here ----
	# ... didn't do that yet, it's mostly irrelevant, except for 
	# onsite variations, e.g. disorder or SSH-like behavior

	# set neighbors (useful for visualization and analysis)
	tb.neighbors = Vector{Vector{Int}}(n)
	for i = 1:n
		tb.neighbors[i] = find(H[i,:])
	end

	return H
end

#= --------------------------------------------------------- =#
function buildHoneycombRestriction{T<:Real}(restrictfun::Function, 
                                   maxscale::T)
	@assert(maxscale>0,"maxscale must be a positive value")

	intralat = [[sqrt(3)*0.5,0.5],[0.0,1.0]]
	shiftmag = 0.5/sqrt(3)
	steplist = -round(1.25*maxscale):round(1.25*maxscale)

	posa = Array{Float64}(0,2) # initialize arrays
	posb = copy(posa)

	# grow A and B sublattice arrays if points inside restriction			
	for n = steplist
		for m = steplist
            newposa = sum(intralat.*[n,m]) + [shiftmag,0]
            newposb = sum(intralat.*[n,m]) - [shiftmag,0]

			if restrictfun(newposa)
				posa = vcat(posa, newposa') 
			end
			if restrictfun(newposb)
				posb = vcat(posb, newposb') 
			end

		end
	end

	# remove dangling bonds
	dangleremoved = 1
	while dangleremoved != 0
		dangleremoved = 0

		# fix A lattice
		delrows=Vector{Int}(0)
		for i in 1:size(posa,1)
			nneighbors = sum(map((x,y) -> interDistance(posa[i,:], [x,y]), posb[:,1], posb[:,2]) .< 0.58)
			if nneighbors < 2 # dangling bond
			    push!(delrows,i)
			end
		end
		posa = posa[setdiff(1:end,delrows),:];
		dangleremoved += length(delrows) 

		# fix B lattice
		delrows=Vector{Int}(0)
		for i in 1:size(posb,1)
			nneighbors = sum(map((x,y) -> interDistance([x,y], posb[i,:]), posa[:,1], posa[:,2]) .< 0.58) 
			if nneighbors < 2 # dangling bond
			    push!(delrows,i)
			end
		end
		posb = posb[setdiff(1:end,delrows),:];
		dangleremoved += length(delrows) 
	end
	
	return vcat(posa, posb)
end


#= --------------------------------------------------------- =#
function solve(tb::TBmodel)
	H = setHopping!(tb)
	ε = eigvals(full(H))

	plot(ε,".")
	return ε
end


#= --------------------------------------------------------- =#
function hopHoneycomb(r1,r2)
	a = interDistance(r1,r2,tb.periodicity)
	t = a>0 && a<.58 ? -1 : 0
end

function hopExponential(r1,r2)
	a = interDistance(r1,r2,tb.periodicity)
	t = a < .15 ? exp(-a/.05) : 0.0
end

f(r,a) = sqrt(sum(r.^2)) < a
f(r) = f(r,15)
