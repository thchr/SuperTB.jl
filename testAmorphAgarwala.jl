include("tb.jl")
close("all")

# set up a model (amorphous system with AII-class hopping)
tb = TBmodel()
tb.pos = rand(576,2)
tb.periodicity=trues(2) # periodic boundary conditions
tb.basisdim = 2
tb.hopfun = hopAgarwala

# solve the system
ε = solve(tb)

# compare with "read-off" data from Agarwala PRL paper (Fig. 2a)
fig2adata = readcsv("fig2a_agarwala.txt")
plot(fig2adata[:,1],fig2adata[:,2],".r")
ylim(-4.5,4.5)

# plot the amorphous lattice (and its neighbors, if for-loops are uncommted)
figure()
#= for site in eachindex(tb.neighbors) # takes a lot of time to plot all 
    for neighbor in tb.neighbors[site] # neighbor connections, unfortunately
        plot(tb.pos[[site,neighbor],1],tb.pos[[site,neighbor],2],"-",color="gray")
    end
end =#
plot(tb.pos[:,1],tb.pos[:,2],".k")
axis("equal")

# calculate the density of states (dos) and plot it
dos(eeig,espan) = map( (x)-> sum(imag(1./(x-eeig-.01im))),espan)
εˢ = linspace(minimum(ε), maximum(ε), 250)
dosˢ = dos(ε, εˢ)

figure()
plot(εˢ, dosˢ, "-k")

