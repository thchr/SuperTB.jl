using Iterators
using ProgressMeter

# we assume that kspan is created by a command like
#    [(i,j) for i=linspace(kxmin,kxmax,Nkx), j=linspace(kymin,kymax,Nky)]
# i.e. that kspan is supplied as an Nkx×Nky array of 2-element tuples
function sdos{T<:Number}(eigval::AbstractArray{S} where S<:Real,
                         eigvec::AbstractArray{T},
                         tb::TBmodel,
                         kspan::AbstractArray=[(i,j) for i=linspace(-3*pi,3*pi,50),j=linspace(-3*pi,3*pi,50)],
                         enespan::Union{Range,AbstractArray}=linspace(minimum(eigval),maximum(eigval),100),
                         eta::Real=0.1
                        )

    eigval_c = eigval .- eta*im # add finite loss

    eigvec_k_abs2 = zeros(Float64,size(eigvec,1))
    expfac = zeros(Complex{Float64},size(tb.pos,1))
    term = zeros(Float64,size(eigval))
    sdosvals = zeros(Float64,(size(kspan)..., size(enespan)...))
    
    @showprogress 2 "Calculating sdos ... " for kindex = CartesianRange(size(kspan))
        expfac .= exp.(im* ( kspan[kindex][1].*tb.pos[:,1] .+
                             kspan[kindex][2].*tb.pos[:,2] ) ) 

        eigvec_k_abs2 .= abs.(squeeze(sum( repmat(expfac,2) .* eigvec, 1),1)).^2

        for eneindex = CartesianRange(size(enespan))
            
            term .= imag( eigvec_k_abs2 ./ (enespan[eneindex] .- eigval_c) ) 
            sdosvals[kindex,eneindex] = - sum(term)
        end
    end

    return sdosvals
end


