type TBmodel
    pos::Array{Real}                    #Position of elements in unit cell
    hop::Array{Array{Complex,1},1}      #Hopping amplitudes to neighbors
    neighbors::Array{Array{Int64,1},1}  #Neighboring indices
    lat::Array{Real}                    #Lattice vectors
end

function interDistance(r::Array,t)

function setHopping!(tb::TBmodel, f::Function, dist::Real)

    

