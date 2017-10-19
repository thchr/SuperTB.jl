
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


tb.pos = rand(sites,2)
tb.hopfun(ri,rj,per) = hopAgarwala(ri, rj, per,
                                   1.0/6,  # R
                                   M,      # M
                                   0.25,   # t₂
                                   0.5,    # λ
                                   1.0/24) # a
# solve the system
H = Hermitian(full(setHopping!(tb)))
eigdecom = eigfact(H)

enespan = linspace(minimum(eigdecom[:values]),maximum(eigdecom[:values]),150)
kx = linspace(-20*pi,20*pi,101); ky = kx
kspan = [(i,j) for i=kx, j=ky]
# compute Bott indices
sdosvals = sdos(eigdecom[:values],
                eigdecom[:vectors],
                tb,
                kspan,
                enespan)



@show size(sdosvals)

figure()
pcolor(kx, enespan, sdosvals[:,51,:].',
       cmap=get_cmap("inferno"),
       vmin=0, vmax=maximum(sdosvals))#surf(enespan,βᵉ)
colorbar()

