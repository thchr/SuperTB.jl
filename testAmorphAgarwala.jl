include("tb.jl")
include("bott.jl")
close("all")

# set up a model (amorphous system with AII-class hopping)
tb = TBmodel()
tb.pos = rand(24^2,2) #rand(576,2)
tb.periodicity=trues(2) # periodic boundary conditions
tb.basisdim = 2
tb.hopfun = hopAgarwala

# solve the system
ε,eigvec = solve(tb; eigvecflag=true)

# plot the eigenenergies
figure()
plot(ε,".")
# compare with "read-off" data from Agarwala PRL paper (Fig. 2a)
fig2adata = readcsv("fig2a_agarwala.txt")
plot(fig2adata[:,1],fig2adata[:,2],".r")
ylim(-4.5,4.5)

# plot the amorphous lattice (and its neighbors, if for-loops are uncommted)
plot(tb; plothopping=false)

# calculate the density of states (dos) and plot it
dos(eeig,espan) = map( (x)-> sum(imag(1./(x-eeig-.01im))),espan)
εˢ = linspace(minimum(ε), maximum(ε), 250)
dosˢ = dos(ε, εˢ)

figure()
plot(εˢ, dosˢ, "-k")


# compute Bott indices
#β = bott(eigvec,tb,size(tb.pos,1))i
bottatindices = 24^2+(-150:150)
β = bott(eigvec,tb,bottatindices)
figure()
plot(ε[bottatindices],β,".-")

