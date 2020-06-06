
include("tb.jl")
include("bott.jl")
include("sdos.jl")
close("all")

# set up a model (amorphous system with AII-class hopping)
tb = TBmodel()
tb.periodicity=trues(2) # periodic boundary conditions
tb.basisdim = 2

# run over this parameter range 
M = -3; 
sites = 24^2 
Nensembles = 20


tb.hopfun(ri,rj,per) = hopAgarwala(ri, rj, per,
                                   1.0/6,  # R
                                   M,      # M
                                   0.25,   # t₂
                                   0.5,    # λ
                                   1.0/24) # a

enespan = linspace(-5,5,150)
kx = linspace(-60*pi,60*pi,301); ky = kx
kspan = [(i,j) for i=kx, j=ky]
sdosvals = zeros((size(kspan)...,size(enespan)...))

for i = 1:Nensembles
    tb.pos = rand(sites,2)
    # solve the system
    H = Hermitian(full(setHopping!(tb)))
    eigdecom = eigfact(H)

    # compute Bott indices
    sdosvals += sdos(eigdecom[:values],
                     eigdecom[:vectors],
                     tb,
                     kspan,
                     enespan)
end
sdosvals /= Nensembles # normalize


figure()
pcolor(kx, enespan, sdosvals[:,151,:].',
       cmap=get_cmap("inferno"),
       vmin=0, vmax=maximum(sdosvals),
       linewidth=0)
colorbar()
PyPlot.svg(true)
savefig("sdos_ensembleaveraged")
