
using SuperTB, PyPlot, LinearAlgebra
close("all")

# set up a model (amorphous system with AII-class hopping)
tb = TBModel(hopAgarwala)
tb.periodicity = trues(2) # periodic boundary conditions
tb.basisdim    = 2

# run over this parameter range 
sites = 24^2 
Nensembles = 20

# spans/preallocation
enespan = range(-5,5,length=150)
kx = range(-60*pi,60*pi,length=301); ky = kx
kspan = [(i,j) for i in kx, j in ky]
sdosvals = zeros(Float64, (size(kspan)...,size(enespan)...))

for i = 1:Nensembles
    tb.pos = rand(sites,2)
    # solve the system
    H = Hermitian(Matrix(setHopping!(tb)))
    F = eigen(H)

    # compute sdos
    sdosvals .+= sdos(F.values, F.vectors, tb, kspan, enespan)
end
sdosvals ./= Nensembles # normalize

figure()
pcolormesh(kx, enespan, sdosvals[:,151,:].',
       cmap=get_cmap("inferno"), vmin=0, vmax=maximum(sdosvals), linewidth=0)
colorbar()
#PyPlot.svg(true)
#savefig("sdos_ensembleaveraged")
