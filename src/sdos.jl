# we assume that kspan is created by a command like
#    [(i,j) for i=range(kxmin,kxmax,lengtH=Nkx), j=range(kymin,kymax,length=Nky)]
# i.e. that kspan is supplied as an Nkx×Nky array of 2-element tuples
function sdos(eigvals::AbstractVector{<:Real},
              eigvecs::AbstractMatrix{<:Number},
              tb::TBModel,
              kspan::AbstractArray = (ks=range(-3π,3π,length=50); [(kx,ky) for kx in ks, ky in ks]),
              enespan::AbstractVector{<:Real} = range(extrema(eigvals)...,length=100),
              eta::Real = 0.1)

    eigvals_c = eigvals .- eta*im # add finite loss

    eigvec_k_abs2 = zeros(Float64,    size(eigvecs, 1))
    expfac        = zeros(ComplexF64, size(tb.pos, 1))
    term          = zeros(Float64,    size(eigvals))
    sdosvals      = zeros(Float64,    (size(kspan)..., size(enespan)...))
 
    @showprogress 2 "Calculating sdos ... " for kidx in CartesianIndices(kspan)
        expfac .= cis.( kspan[kidx][1].*tb.pos[:,1] .+ kspan[kidx][2].*tb.pos[:,2] ) 

        eigvec_k_abs2 .= abs2.(dropdims(sum(repeat(expfac, tb.basisdim) .* eigvecs, dims=1), 
                                        dims=1)) # TODO: ... probably allocates...

        for eneidx = eachindex(enespan)
            term .= imag.( eigvec_k_abs2 ./ (enespan[eneidx] .- eigvals_c) ) 
            sdosvals[kidx,eneidx] = -sum(term)
        end
    end

    return sdosvals
end


