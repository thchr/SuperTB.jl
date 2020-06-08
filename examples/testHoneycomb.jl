using SuperTB, PyPlot

plt.close("all")

tb = TBModel(hopHoneycomb)
tb.pos = buildHoneycombRestriction(x->circularRestriction(x,10), 10)
tb.periodicity = falses(2)

plot(tb, plothopping=true)

# crude dos
ε = solve(tb)
dos(eeig,espan) = map(e -> sum(imag(inv.((e - .01im) .- eeig))), espan)
εˢ = range(minimum(ε), maximum(ε), length=250)

figure()
plot(εˢ, dos(ε, εˢ), "-k")