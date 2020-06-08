using SuperTB, PyPlot, ProgressMeter, LinearAlgebra
close("all")

# set up a model (amorphous system with AII-class hopping)
tb = TBModel(hopAgarwala)
tb.periodicity = trues(2) # periodic boundary conditions
tb.basisdim    = 2

# run over this parameter range 
M = range(-4.5, 4.5, length=50)
sites = 12:12:24^2 
Nsamples = 1 # number of samples in ensemble

β = zeros(length(M), length(sites), Nsamples)

@showprogress 2 "Building phase diagram (#sites-progress)..." for s = eachindex(sites)
    for g = 1:Nsamples
        tb.pos = rand(sites[s],2)
   
        for m = eachindex(M)
            tb.hopfun(ri,rj,per) = hopAgarwala(ri, rj, per,
                                               1.0/6,  # R
                                               M[m],   # M
                                               0.25,   # t₂
                                               0.5,    # λ
                                               1.0/24) # a
            # solve the system
            H = Hermitian(Matrix(setHopping!(tb)))
            F = eigen(H)

            # compute Bott indices
            β[m,s,g] = bott(F.vectors,tb,sites[s])[1]
        end
    end
end
βᵉ = dropdims(sum(β,3),dims=3)./Nsamples
ρ = sites./24^2

figure()
surf(ρ,M,βᵉ)
