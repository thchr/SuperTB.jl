include("tb.jl")
close("all")

tb = TBmodel()
tb.pos = buildHoneycombRestriction(x->f(x,10), 10)
tb.periodicity=falses(2)
tb.hopfun = hopHoneycomb

figure()
plot(tb.pos[:,1],tb.pos[:,2],".k")
axis("equal")

ε = solve(tb)

dos(eeig,espan) = map( (x)-> sum(imag(1./(x-eeig-.01im))),espan)
εˢ = linspace(minimum(ε), maximum(ε), 250)

dosˢ = dos(ε, εˢ)

figure()
plot(εˢ, dosˢ, "-k")

