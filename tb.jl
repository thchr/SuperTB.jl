import SparseSuite

type TBmodel
    pos::Vector{Vector{<:Real}}           # position of elements in unit cell
    hop::Vector{Vector{<:Complex}}        # hopping amplitudes to neighbors
    neighbors::Vector{Vector{<:Integer}}  # neighboring indices
    lat::Vector{Vector{<:Real}}           # lattice vectors
end

function interDistance(pos::Vector{Vector{<:Real}}; isperiodic::Vector{Int}=[0,0,0])
	n = length(pos)
	D = zeros(eltype(pos[1]),n,n)

	if iszero(isperiodic) # not periodic in any dimension
		for i in eachindex(pos)
			for j = (i+1):n # fill upper triangle with inter-distances
				D[i,j] = sqrt( sum( (pos[i] - pos[j]).^2 ) )
			end
		end

	else # periodic in some directions (specifically, distance on an n-torus, with n = sum(isperiodic))
		xyzdims = length(pos[1])
		for i in eachindex(pos)
			for j = (i+1):n # fill upper triangle with inter-distances
				s = 0.0;
				rij = abs.(pos[i]-pos[j])
				for xyz = 1:xyzdims
					if isperiodic[xyz] # dimension is periodic
						s += min(rij[xyz],1-rij[xyz])^2
					else			   # dimension is not periodic
						s += rij[xyz]^2
                    end
					D[i,j] = sqrt(s)
				end
			end
		end
	end

	return Symmetric(D, :U)
end

function setHopping!(tb::TBmodel, f::Function, dist::Real)

    
end
