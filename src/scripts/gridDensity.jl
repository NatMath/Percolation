# What is the effect of density on the number of sites adjacent to things
# Expects size as dim tupple, so if 2D, want (L,L)
function selfAdjacency(size, repeats = 10)
    ρs = 0:0.01:1
    counts = similar(ρs)
    for i in 1:length(ρs)
        ρ = ρs[i]
        count = round(Int, prod(size) * ρ)
        vect = zeros(prod(size))
        vect[1:count] .= 1.
        countBuffer = 0.
        for k in 1:repeats
            grid = reshape(vect, size)
            shuffle!(grid)

            adjacents = findNeighboured(grid, 1.)
            countBuffer += length(adjacents)
        end
        counts[i] = countBuffer / repeats
    end
    println("Plotting")
    plot(ρs, counts, xlabel = "Site Density", ylabel = "Number of Sites with Adjacency")
end
export selfAdjacency
