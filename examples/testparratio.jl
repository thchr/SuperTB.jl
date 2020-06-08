using SuperTB, PyPlot, ProgressMeter, LinearAlgebra

# set up a model (amorphous system with AII-class hopping)
hopfun(ri,rj,per) = hopAgarwala(ri, rj, per,
                                1.0/6,  # R
                                -.5,    # M
                                0.25,   # t₂
                                0.5,    # λ
                                1.0/64) # a
tb = TBModel(hopfun)
tb.periodicity = trues(2) # periodic boundary conditions
tb.basisdim    = 2

# run over this parameter range 
sitessqrt  = 24
sites      = sitessqrt^2 
Nensembles = 5
getbottatindex = sites .+ (-(2*sitessqrt+1):(2*sitessqrt+1))

# preallocation
par  = zeros(sites*tb.basisdim, Nensembles)      # participation ratios
εeig = zeros(sites*tb.basisdim, Nensembles)      # energy eigenvalues
β    = zeros(size(getbottatindex,1), Nensembles) # bott indices

close("all")
# compute and visualize during...
@showprogress 1 "Calculating ensembles..." for i in 1:Nensembles
    tb.pos = rand(sites,2)
    # solve the system
    H = Hermitian(Matrix(setHopping!(tb)))
    F = eigen(H)

    # compute participation ratios & plot
    par[:,i]  = parratio(F.vectors)
	εeig[:,i] = F.values
	figure()
	plot(εeig[:,i], par[:,i], ".")
        
    # compute bott indices & plot
    β[:,i] = bott(F.vectors,  tb,getbottatindex)
	plot(εeig[getbottatindex,i], β[:,i]./10, ".") # scale to fit in same plot...
end

#PyPlot.svg(true)
#savefig("sdos_ensembleaveraged")