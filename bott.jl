using Base.LinAlg
using ProgressMeter

function bott{T<:Number}(eigvec::AbstractArray{T},
                         tb::TBmodel,
                         outputRange::Union{UnitRange,Integer}=1:size(eigvec,2)
                         )

    θx = tb.pos[:,1].*(2*pi) # assuming x and y to be between 0 and 1
    θy = tb.pos[:,2].*(2*pi)

    expx = exp.(im*θx)
    expy = exp.(im*θy)
    expx = repmat(expx, tb.basisdim)
    expy = repmat(expy, tb.basisdim)

    if size(eigvec,2) > maximum(outputRange)
        eigvec = view(eigvec,:,1:maximum(outputRange))
    end

    X = eigvec'*(expx.*eigvec)
    Y = eigvec'*(expy.*eigvec)

    β = zeros(Float64, size(outputRange))
    @showprogress 2 "Computing Bott indices ..." for n = outputRange
        x = view(X,1:n,1:n)
        y = view(Y,1:n,1:n)
        β[n-outputRange[1]+1] = imag(sum(log.(eigvals(y*x*y'*x'))))/(2*pi)
    end

    return β
end
