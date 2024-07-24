# Counts the number of direct neighbours of each cells are the target value and that match predicate
function countNeighbours(grid, target, predicate = x -> true) # for 2D grid
    neighbours = fill!(similar(grid), 0.)
    neighbours[begin:end-1, :] .+= (grid[begin+1:end, :] .== target)
    neighbours[begin+1:end, :] .+= (grid[begin:end-1, :] .== target)
    neighbours[:, begin:end-1] .+= (grid[:, begin+1:end] .== target)
    neighbours[:, begin+1:end] .+= (grid[:, begin:end-1] .== target)
    neighbours[findall(x -> !predicate(x), grid)] .= 0.
    return neighbours
end
export countNeighbours


# Returns 1D list of indices
# Returns all the cell indices where there are neighbours of type target, and the position value matches predicate
findNeighboured(grid, target, predicate = x -> true) = findall(x -> x>0, countNeighbours(grid, target, predicate))
export findNeighboured

function spatialCorrelation(grid, cutoffDist = size(grid, 1) - 1)
    offsets = [(x,y) for x in 0:cutoffDist, y in 0:cutoffDist if x^2+y^2 <= cutoffDist^2]
    # Set of vectors which contain the behaviour as function of d. Indexed by distance^2+1.
    # For example, offset (3,4) gives index 3^2+4^2+1 = 26. Offset (5,0) also gives index 26, because they are at the same distance
    cummAB = zeros(cutoffDist^2 + 1)
    cummA =  zeros(cutoffDist^2 + 1)
    cummB =  zeros(cutoffDist^2 + 1)
    counts = zeros(Int, cutoffDist^2 + 1)
    buffer = similar(grid)
    for (xoff, yoff) in offsets
        A = @view grid[begin+xoff:end, begin+yoff:end]
        B = @view grid[begin:end-xoff, begin:end-yoff]
        bufferView = @view buffer[begin:end-xoff, begin:end-yoff]
        # All the points of A,B are currently at a distance sqrt(d^2) from each other. 
        ind = xoff^2 + yoff^2 + 1 #Index

        bufferView .= A .* B
        cummAB[ind] += sum(bufferView)
        cummA[ind]  += sum(A)
        cummB[ind]  += sum(B)
        counts[ind] += prod(size(A))
    end
    autocorrs = @. cummAB/counts - cummA*cummB/counts^2
    distances = sqrt.(0:cutoffDist^2)
    return (distances[counts .!= 0], autocorrs[counts .!= 0])
end
export spatialCorrelation