function plot(tb::TBModel; plothopping::Bool=false)
	figure()
	if plothopping
		x_cncts = Vector{Float64}[]
		y_cncts = similar(x_cncts)
		for site in eachindex(tb.neighbors)
			for neighbor in tb.neighbors[site]
				push!(x_cncts, tb.pos[[site,neighbor],1])
				push!(y_cncts, tb.pos[[site,neighbor],2])
			end
		end
		# faster to collect all the connections in one array than to issue multiple calls to 
		# PyPlot for each connection
		plot(hcat(x_cncts...), hcat(y_cncts...), "-", color="gray")
		# TODO: Avoid plotting same connection twice...
	end

	plot(tb.pos[:,1], tb.pos[:,2], ".k")
	axis("equal")
end