
using SuperTB, PyPlot, DelimitedFiles
close("all")

# set up a model (amorphous system with AII-class hopping)
hop = (ri, rj, per) -> hopAgarwala(ri, rj, per)
tb = TBModel(hop)
tb.pos = rand(24^2,2) #rand(576,2)
tb.periodicity=trues(2) # periodic boundary conditions
tb.basisdim = 2

# solve the system
ε, eigvec = solve(tb; eigvecflag=true)

# plot the eigenenergies
figure()
plot(ε, ".")
# compare with manually read-off data from Agarwala PRL paper (Fig. 2a)
fig2adata = readdlm(join(pkgdir(SuperTB)*"/data/fig2a_agarwala.txt"), ',', Float64)
plot(fig2adata[:,1], fig2adata[:,2], ".r")
ylim(-4.5, 4.5)

# plot the amorphous lattice (and its neighbors, if for-loops are uncommted)
plot(tb)

# calculate the density of states (dos) and plot it
dos(eeig, espan) = map( (x)-> sum(imag(inv.((x - 0.01im) .- eeig))), espan)
εˢ = range(minimum(ε), maximum(ε), length=250)
dosˢ = dos(ε, εˢ)

figure()
plot(εˢ, dosˢ, "-k")

# compute Bott indices
bottatindices = size(tb.pos,1) .+ (-150:150) # focus on the middle part of the spectrum...
β = bott(eigvec, tb, bottatindices)
figure()
plot(ε[bottatindices],β,".-")