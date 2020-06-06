include("tb.jl")
include("bott.jl")
include("parratio.jl")
close("all")

# set up a model (amorphous system with AII-class hopping)
tb = TBmodel()
tb.periodicity=trues(2) # periodic boundary conditions
tb.basisdim = 2

# run over this parameter range 
M = -.5; 
sitessqrt = 48; 
sites = sitessqrt^2 
Nensembles = 5


tb.hopfun(ri,rj,per) = hopAgarwala(ri, rj, per,
                                   1.0/6,  # R
                                   M,      # M
                                   0.25,   # t₂
                                   0.5,    # λ
                                   1.0/64) # a

getbottatindex=(sites-(2*sitessqrt+1)):(sites+(2*sitessqrt+1))
par = zeros(sites*tb.basisdim,Nensembles)
eneeig = zeros(sites*tb.basisdim,Nensembles)
β = zeros(size(getbottatindex,1),Nensembles)
@showprogress 1 "Calculating ensembles..." for i = 1:Nensembles
    tb.pos = rand(sites,2)
    # solve the system
    H = Hermitian(full(setHopping!(tb)))
    eigdecom = eigfact(H)

    # compute Bott indices
    par[:,i] = parratio(eigdecom[:vectors])
	eneeig[:,i] = eigdecom[:values]

	
	figure()
	plot(eneeig[:,i],par[:,i],".")
		
	β[:,i] = bott(eigdecom[:vectors],tb,getbottatindex)
	plot(eneeig[getbottatindex,i],β[:,i]/10,".")
end





#PyPlot.svg(true)
#savefig("sdos_ensembleaveraged")
