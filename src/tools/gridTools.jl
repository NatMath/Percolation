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

function findClusters2D(grid::Matrix{Bool})
    # Based on the naive Hoshen-Kopelman algorithm
    #https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html

    # Expects a grid of booleans values, with true being occupied
    # Doesn't need to be square

    L = size(grid, 1)
    N = prod(size(grid))

    labels = collect(1:N)

    for i = 1:size(grid, 2), j = 1:size(grid, 1)
        idx = (i-1)*L+j
        grid[j,i] || continue #If unoccupied, continue to next grid cell
        left = (i > 1 && grid[j, i-1]) 
        up = (j > 1 && grid[j-1, i])   

        !left && !up && continue
        if left && ! up
            labels[idx] = labels[idx-L]
        elseif !left && up
            labels[idx] = labels[idx-1]
        elseif left && up
            union!(labels, idx-L, idx-1)
            labels[idx] = labels[idx-L]
        end
    end
    return [find!(labels, x) for x in 1:N]
end
export findClusters2D
# Helper methods. Not exported
function union!(labels, x, y)
    labels[find!(labels, x)] = find!(labels, y)
    return labels
end
function find!(labels, x)
    y = x
    while labels[y] != y #Find the root
        y = labels[y]
    end 
    while labels[x] != x #Set all the nodes on the way directly to the root
        z = labels[x]
        labels[x] = y
        x = z
    end
    return y
end

function countClusters2D(grid::Matrix{Bool})
    clusters = reshape(findClusters2D(grid), size(grid)) .* grid
    return length(Set(clusters)) - 1 #Remove 0, which corresponds to the points where the grid is unoccupied
end
export countClusters2D

function clusterSizes2D(grid::Matrix{Bool})
    clusters = reshape(findClusters2D(grid), size(grid)) .* grid
    sizes = [count(==(n), clusters) for n in unique(clusters[clusters .!= 0])]
    return sizes
end
export clusterSizes2D