#= --------------------------------------------------------- =#
mutable struct TBModel{T<:Real}
    pos::Matrix{T}  		# position of elements in unit cell
    periodicity::BitVector  # periodicity along latx,laty,latz (boolean = yes/no)
    neighbors::Vector{Vector{Int}} # neighboring indices
	hopfun::Function        # hopping function (input is two positions; output is basisdim×basisdim array)
    basisdim::Int           # dimension of basis on each site (must be same everywhere, for now)
end
TBModel(hopfun) = TBModel( Matrix{Float64}(undef,0,2),  # pos
                     falses(2),                   # periodicity
					 Vector{Vector{Int}}(),       # neighbors
					 hopfun,                      # hopfun
                     1)                           # basisdim


#= --------------------------------------------------------- =#
function interDistance(ri::AbstractVector{<:Real}, rj::AbstractVector{<:Real},
					   periodicity::AbstractVector{Bool}=falses(length(first(ri))))

	# TODO: get rid of all the array allocations here...!
	# not periodic in any dimension; regular flat metric
	if all(!, periodicity)
		return sqrt( sum( abs2.(ri - rj) ) )

	# periodic in some directions; use distance on a flat n-torus,
	# with n = sum(periodicity)
	else 
		# TODO: don't hardcode this to a [0,1]² box
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
function interAngle(ri::AbstractVector{<:Real}, rj::AbstractVector{<:Real},
                    periodicity::AbstractVector{Bool}=falses(length(first(ri))))
    #= calculates the angle between two 2D points ri and rj, on a
       possibly periodic (flat) manifold =#
    rij = ri - rj
    for xy = 1:2 # fold to shortest distance over periodic directions
		if periodicity[xy] 
			# TODO: Don't hardcode this to a [0,1]² box...
            if abs(rij[xy]) > abs(rij[xy]-1)
                rij[xy] -= 1
            elseif abs(rij[xy]) > abs(rij[xy]+1)
                rij[xy] += 1
            end
        end
    end
    return atan(rij[2],rij[1])
end


#= --------------------------------------------------------- =#
function setHopping!(tb::TBModel)

	n = size(tb.pos,1)
    basisiter = Base.OneTo(tb.basisdim)
    H = spzeros(ComplexF64, n*tb.basisdim, n*tb.basisdim)

	#= fill upper triangle with hopping amplitude by evaluating
    the function hopfun as a function of (i,j)-positions  =#
	for i = Base.OneTo(n)
		for j = (i+1):n
			t = tb.hopfun(tb.pos[i,:], tb.pos[j,:], tb.periodicity)
            for k = basisiter
                for l = basisiter
                    if !iszero(t[k,l])
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
	tb.neighbors = Vector{Vector{Int}}(undef, n)
	for i = 1:n
		# TODO: This is a horrible way to do this with a sparse array, lol... Should use a 
		#		combination of column iteration and nzrange instead...
		tb.neighbors[i] = findall(!iszero, H[i,1:n])
	end

	return H
end

#= --------------------------------------------------------- =#
function buildHoneycombRestriction(restrictfun::Function, maxscale::Real)
	maxscale≤0 && throw(DomainError(maxscale, "maxscale must be a positive value"))

	intralat = [ [sqrt(3)*0.5, 0.5], [0.0, 1.0] ]
	shiftmag = 0.5/sqrt(3)
	steplist = -round(1.25*maxscale):round(1.25*maxscale)

	posa = Matrix{Float64}(undef, 0,2) # initialize arrays
	posb = similar(posa)

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
		delrows = Int[]
		for i in 1:size(posa,1)
			nneighbors = sum(map((x,y) -> interDistance(posa[i,:], [x,y]), posb[:,1], posb[:,2]) .< 0.58)
			if nneighbors < 2 # dangling bond
			    push!(delrows,i)
			end
		end
		posa = posa[setdiff(1:end,delrows),:];
		dangleremoved += length(delrows) 

		# fix B lattice
		delrows = Int[]
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
function solve(tb::TBModel; eigvecflag::Bool=false)
	H = setHopping!(tb)
    if !eigvecflag
    	ε = eigvals(Hermitian(Matrix(H)))
        return ε
    else
        F = eigen(Hermitian(Matrix(H)))
        return F.values, F.vectors
    end
	
end


#= --------------------------------------------------------- =#
function hopHoneycomb(ri::AbstractVector{<:Real}, rj::AbstractVector{<:Real},
                      periodicity::AbstractVector{Bool}=falses(length(first(ri))))
    #= Nearest neighbor hopping on a honeycomb (graphene-like) lattice 
       with unity hopping amplitude =#
	a = interDistance(ri, rj, periodicity)
	t = a<.58 ? -1.0 : 0.0

end

#= --------------------------------------------------------- =#
function hopAgarwala(ri::AbstractVector{<:Real}, rj::AbstractVector{<:Real},
                     periodicity::AbstractVector{Bool}=trues(length(first(ri))),
                     R=1.0/6, M=-0.5, t₂=0.25, λ=0.5, a=1.0/24)
    #= Implementation of the hopping type described in Eq. (5) of Agarwala & 
       Shenoy, PRL 118, 236402 (2017) =#

	r = interDistance(ri, rj, periodicity)
    
    if r > R
		t = zeros(ComplexF64, 2, 2)
    elseif r > zero(eltype(r))
        ϕ = interAngle(ri, rj, periodicity) # works properly also in periodic cases

        # the minus sign in front of exp _must_ be there in front of _both_
        # offdiag terms in order to get the right result; i.e. t is not Hermitian 
		# on its own (but H is still Hermitian, by construction)
		imcisϕ = -cis(-ϕ)*im
		sin²ϕ  = sin(ϕ)^2
        offdiagᴬᴮ = ( imcisϕ + λ*(sin²ϕ*(1+im)-1) ) / 2
        offdiagᴮᴬ = ( imcisϕ + λ*(sin²ϕ*(1-im)-1) ) / 2 
        diagᴬᴬ = 0.5 * (-1 + t₂)
        diagᴮᴮ = 0.5 * (+1 + t₂)

        t = exp(1 - r/a) .* [ diagᴬᴬ offdiagᴬᴮ; offdiagᴮᴬ diagᴮᴮ ]
    else # on-site/self term (r=0)
        t = [(2+M) (1-im)*λ; (1+im)*λ -(2+M)] # this term must be Hermitian though
    end

	return t
end

#= --------------------------------------------------------- =#
circularRestriction(r,a) = sqrt(sum(r.^2)) < a
circularRestriction(r) = circularRestriction(r,15)
