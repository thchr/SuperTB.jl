using Base.LinAlg

function parratio{T<:Number}(eigvec::AbstractArray{T})

	denom = sum(abs.(eigvec).^2,2).^2
	numer = sum(abs.(eigvec).^4,2)
	return numer./denom

end
