module SuperTB

using LinearAlgebra
using ProgressMeter
using SparseArrays

include("tb.jl")
export TBModel, setHopping!, solve,
    buildHoneycombRestriction, hopHoneycomb, 
    hopAgarwala, circularRestriction

include("bott.jl")
export bott

include("sdos.jl")
export sdos

include("participation.jl")
export parratio

import PyPlot: plot
using PyPlot
include("plot.jl")
export plot

end # module