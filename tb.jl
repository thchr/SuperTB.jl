using PyPlot
import Base.Math

#= --------------------------------------------------------- =#
mutable struct TBmodel{T<:Real}
    pos::AbstractArray{T,2}  # position of elements in unit cell
    periodicity::AbstractVector{Bool} # periodicity along latx,laty,latz (boolean = yes/no)
    neighbors::AbstractArray{Q} where Q <: AbstractArray{S} where S<: Integer # neighboring indices
	hopfun::Function 		 # hopping function (input is two positions; output is basisdim×basisdim array)
    basisdim::Integer        # dimension of basis on each site (must be same everywhere, for now)
#    lat::AbstractVector{AbstractVector{T}}       # lattice vectors
end
TBmodel() = TBmodel( Array{Float64}(0,2),      # pos
                     falses(2),                # periodicity
					 Vector{Vector{Int}}(0),   # neighbors
					 (x,y) -> zero(eltype(x)), # hopfun
                     1 )                       # basisdim


#= --------------------------------------------------------- =#
function interDistance{T<:Real}(ri::AbstractVector{T},
                                rj::AbstractVector{T},
                                periodicity::AbstractVector{Bool}=[false,false])

	# not periodic in any dimension; regular flat metric
	if iszero(periodicity) 
		return sqrt( sum( (ri - rj).^2 ) )

	# periodic in some directions; use distance on a flat n-torus,
	# with n = sum(periodicity)
	else 
		s = zero(eltype(ri))
		rij = abs.(ri - rj)
		for xyz in eachindex(ri)
			if periodicity[xyz] # dimension is periodic
				s += min(rij[xyz], 1-rij[xyz])^2
			else			    # dimension is not periodic
				s += rij[xyz]^2
			end
		end
		return sqrt(s)
	end

end


#= --------------------------------------------------------- =#
function interAngle{T<:Real}(ri::AbstractVector{T},
                             rj::AbstractVector{T},
                             periodicity::AbstractVector{Bool}=[false,false])
    #= calculates the angle between two 2D points ri and rj, on a
       possibly periodic (flat) manifold =#
    rij = ri - rj
    for xy = 1:2 # fold to shortest distance over periodic directions
        if periodicity[xy] 
            if abs(rij[xy]) > abs(rij[xy]-1)
                rij[xy] -= 1
            elseif abs(rij[xy]) > abs(rij[xy]+1)
                rij[xy] += 1
            end
        end
    end
    return atan2(rij[2],rij[1])
end


#= --------------------------------------------------------- =#
function setHopping!(tb::TBmodel)
	n = size(tb.pos,1)
    basisiter = 1:tb.basisdim
    H = spzeros(Complex{Float64}, n*tb.basisdim, n*tb.basisdim)

	#=  fill upper triangle with hopping amplitude by evaluating
    the function hopfun as a function of (i,j)-positions  =#
	for i = 1:n
		for j = (i+1):n
			t = tb.hopfun(tb.pos[i,:], tb.pos[j,:], tb.periodicity)
            for k = basisiter
                for l = basisiter
                    if t[k,l] != zero(eltype(t))
                        H[i+(k-1)*n,j+(l-1)*n] = t[k,l]  # upper triangular part
                    end
                end
			end
		end
	end
   
    #= fill out the lower triangle, using Hermicity of H =#
    H = H + H' # when done this way, it has be done _before_ diagonal part is entered
    
	# ---- diagonal of H should be computed here ----
    for i = 1:n
        # stupid implementation; but easier to extend to index-dependence
        t = tb.hopfun(tb.pos[i,:], tb.pos[i,:], tb.periodicity) 
        for k = basisiter
            for l = basisiter
                H[i+(k-1)*n,i+(l-1)*n] = t[k,l]
            end
        end
    end

	# set neighbors (useful for visualization and analysis)
	tb.neighbors = Vector{Vector{Int}}(n)
	for i = 1:n
		tb.neighbors[i] = find(H[i,1:n])
	end

	return H
end

#= --------------------------------------------------------- =#
function buildHoneycombRestriction{T<:Real}(restrictfun::Function, 
                                            maxscale::T)
	@assert(maxscale>0,"maxscale must be a positive value")

	intralat = [ [sqrt(3)*0.5, 0.5], [0.0, 1.0] ]
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
function solve(tb::TBmodel; eigvecflag::Bool=0)
	H = setHopping!(tb)
    if !eigvecflag
    	ε = eigvals(full(H))
        return ε
    else
        (ε,eigvec) = eig(full(H))
        return ε, eigvec
    end
	
end


#= --------------------------------------------------------- =#
function hopHoneycomb{T<:Real}(ri::AbstractVector{T}, rj::AbstractVector{T},
                               periodicity::AbstractVector{Bool})
    #= Nearest neighbor hopping on a honeycomb (graphene-like) lattice 
       with unity hopping amplitude =#
	a = interDistance(ri,rj,periodicity)
	t = a>0 && a<.58 ? -1 : 0

end

#= --------------------------------------------------------- =#
function hopAgarwala{T<:Real}(ri::AbstractVector{T}, rj::AbstractVector{T},
                              periodicity::AbstractVector{Bool},
                              R=1.0/6, M=-0.5, t₂=0.25, λ=0.5, a=1.0/24)
    #= Implementation of the hopping type described in Eq. (5) of Agarwala & 
       Shenoy, PRL 118, 236402 (2017) =#

	r = interDistance(ri, rj, periodicity)
    
    if r > R
        t = zeros(Complex{Float64}, 2, 2)
    elseif r > zero(eltype(r))
        ϕ = interAngle(ri, rj, periodicity) # will work properly also in periodic cases

        # the minus sign in front of exp _must_ be there in front of _both_
        # offdiag terms in order to get the right result; i.e. t is not Hermitian 
        # on its own (but H is still Hermitian, by construction)
        offdiagᴬᴮ = ( - exp(-ϕ*im)*im + λ*(sin(ϕ)^2*(1+im)-1) ) / 2
        offdiagᴮᴬ = ( - exp(+ϕ*im)*im + λ*(sin(ϕ)^2*(1-im)-1) ) / 2 
        diagᴬᴬ = 0.5 * (-1 + t₂)
        diagᴮᴮ = 0.5 * (+1 + t₂)

        t = e * exp(-r/a) .* [ diagᴬᴬ offdiagᴬᴮ; offdiagᴮᴬ diagᴮᴮ ]
    else
        t = [(2+M) (1-im)*λ; (1+im)*λ -(2+M)] # this term must be Hermitian though
    end

end

#= --------------------------------------------------------- =#
function plot(tb::TBmodel; plothopping::Bool=0)
	figure()
	if plothopping
		for site in eachindex(tb.neighbors)    # takes a lot of time to plot all 
			for neighbor in tb.neighbors[site] # neighbor connections, unfortunately
				plot(tb.pos[[site,neighbor],1],tb.pos[[site,neighbor],2],"-",color="gray")
			end
		end 
	end
	plot(tb.pos[:,1],tb.pos[:,2],".k")
	axis("equal")
end



#= --------------------------------------------------------- =#
circularRestriction(r,a) = sqrt(sum(r.^2)) < a
circularRestriction(r) = circularRestriction(r,15)
