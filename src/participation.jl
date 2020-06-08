"""
	parratio(eigvec::AbstractMatrix{<:Number})

Compute participation ratio of an array of eigenvectors, supplied as a Matrix
"""
function parratio(eigvec::AbstractMatrix{<:Number})
	denom = sum(abs2, eigvec, dims=2).^2
	numer = sum(x->abs2(x)^2, eigvec,dims=2)

	return numer./denom
end
