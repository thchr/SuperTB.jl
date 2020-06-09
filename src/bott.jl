function bott(eigvec::AbstractMatrix{<:Number}, tb::TBModel,
              outputRange::AbstractRange=Base.OneTo(size(eigvec,2)))

    # computes exp(iθx) & exp(iθy) with θx ≡ x/2π & θy ≡ y/2π and [x y] = tb.pos[:,1:2],
    # and then block-repeats it for cases with non-unity basis sizes (assumed split over a
    # similar block structure). This implicitly assumes x and y to be ∈[0,1]
    expx = repeat(cis.(tb.pos[:,1].*(2π)), tb.basisdim)
    expy = repeat(cis.(tb.pos[:,2].*(2π)), tb.basisdim)

    if outputRange != 1:size(eigvec, 2)
        eigvec = @view eigvec[:,outputRange]
    end

    X = eigvec'*(expx.*eigvec)
    Y = eigvec'*(expy.*eigvec)

    β = zeros(Float64, length(outputRange))
    @showprogress 2 "Computing Bott indices ..." for idx in eachindex(outputRange)
        # Not worth it to do views here, I think; BLAS is more important than GC...
        X′ = X[Base.OneTo(idx), Base.OneTo(idx)]
        Y′ = Y[Base.OneTo(idx), Base.OneTo(idx)]
        β[idx] = imag(sum(log.(eigvals(Y′*X′*Y′'*X′'))))/(2*pi)
    end

    return β
end
