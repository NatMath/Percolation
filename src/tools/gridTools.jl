# Counts the number of direct neighbours of each cells are the target value and that match predicate
# Multi-dimensional and more performant than previous due to views instead of slices
function countNeighbours(grid, target, predicate = Returns(true))
    neighbours = zeros(Int8, size(grid))
    dims = ndims(grid)
    for d in 1:dims
        b = firstindex(grid, d) #Gives access to begin/end indexing outside of brackets
        e = lastindex(grid, d)
        selectdim(neighbours, d, b:(e-1)) .+= (selectdim(grid, d, (b+1):e) .== target)
        selectdim(neighbours, d, (b+1):e) .+= (selectdim(grid, d, b:(e-1)) .== target)
    end
    neighbours[findall(!predicate, grid)] .= 0
    return neighbours
end
export countNeighbours

# Returns 1D list of indices
# Returns all the cell indices where there are neighbours of type target, and the position value matches predicate
findNeighboured(grid, target, predicate = Returns(true)) = findall(x -> x>0, countNeighbours(grid, target, predicate))
export findNeighboured

function findClusters(grid::Matrix{Bool}, ::NormalBoundary)
    # Based on the naive Hoshen-Kopelman algorithm
    # https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html

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
    return reshape([find!(labels, x) for x in 1:N], size(grid))
end

function findClusters(grid::Matrix{Bool}, ::PeriodicBoundary)
    L = size(grid, 1)
    N = prod(size(grid))
    labels = collect(1:N)

    for i = 1:size(grid, 2), j = 1:size(grid, 1)
        # Same as for normal boundary conditions
        idx = (i-1)*L+j
        grid[j,i] || continue #If unoccupied, continue to next grid cell
        left = (i > 1 && grid[j, i-1]) 
        up = (j > 1 && grid[j-1, i])   

        if left && ! up
            labels[idx] = labels[idx-L]
        elseif !left && up
            labels[idx] = labels[idx-1]
        elseif left && up
            union!(labels, idx-L, idx-1)
            labels[idx] = labels[idx-L]
        end
        # On the right side of the grid, and the first cell on the left is occupied
        # => Connecting to a cluster due to periodic boundary
        if i == size(grid, 2) && grid[j, 1]
            union!(labels, idx, j)
        end
        if j == size(grid, 1) && grid[1, i] # Connecting through the bottom
            union!(labels, idx, (i-1)*L+1)
        end
    end
    return reshape([find!(labels, x) for x in 1:N], size(grid))
end
export findClusters

# Helper methods. Not exported
# Set the target of the root of one cluster to the root of the other
# maintains the guarantee that if labels[a]=b, b<=a. Labels "point down" a chain until the root
function union!(labels, x, y) 
    rootx = find!(labels, x)
    rooty = find!(labels, y)
    labels[max(rootx, rooty)] = min(rootx, rooty) #Set the root of the x cluster to point to the root of the y cluster
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

# Multi-dimensional versions of the cluster finding algorithm. 
# Slightly slower than the previous due to the use of union! if there is only one neighbour. Might be improves, but not significantly
# Keep both with dispatch on Matrix vs Array
function findClusters(grid::Array{Bool}, ::NormalBoundary)
    N = prod(size(grid))
    labels = collect(1:N)
    linearIndices = LinearIndices(grid)
    cartesianIndices = CartesianIndices(grid)

    for idx in cartesianIndices
        lidx = linearIndices[idx]
        grid[idx] || continue # If unoccupied, continue to next

        for d in 1:ndims(grid)
            idx[d] == firstindex(grid, d) && continue # Currently at the first index along d
            I = Base.setindex(zero(idx), 1, d) # Cartesian index that is all zero, except 1 along d
            if grid[idx - I]
                union!(labels, lidx, linearIndices[idx-I])
            end
        end
    end
    return reshape([find!(labels, x) for x in 1:N], size(grid))
end
function findClusters(grid::Array{Bool}, ::PeriodicBoundary)
    N = prod(size(grid))
    labels = collect(1:N)
    linearIndices = LinearIndices(grid)
    cartesianIndices = CartesianIndices(grid)

    for idx in cartesianIndices
        lidx = linearIndices[idx]
        grid[idx] || continue # If unoccupied, continue to next

        for d in 1:ndims(grid)
            # Same as normal boundary conditions
            idx[d] == firstindex(grid, d) && continue # Currently at the first index along d
            I = Base.setindex(zero(idx), 1, d) # Cartesian index that is all zero, except 1 along d
            if grid[idx - I]
                union!(labels, lidx, linearIndices[idx-I])
            end
            
            #Periodic boundary condition
            if idx[d] == lastindex(grid, d) && grid[Base.setindex(idx, firstindex(grid, d), d)]
                # Last cell along the dimension and the next cell over is occupied
                union!(labels, lidx, linearIndices[Base.setindex(idx, firstindex(grid, d), d)])
            end
        end
    end
    return reshape([find!(labels, x) for x in 1:N], size(grid))
end

# Works for any boundary conditions
function countClusters(grid::Array{Bool}, bc::BoundaryCondition)
    clusters = findClusters(grid, bc) .* grid
    return length(Set(clusters)) - 1 #Remove 0, which corresponds to the points where the grid is unoccupied
end
export countClusters

function clusterSizes(grid::Array{Bool}, bc::BoundaryCondition)
    clusters = findClusters(grid, bc) .* grid
    sizes = [count(==(n), clusters) for n in unique(clusters[clusters .!= 0])]
    return sizes
end
export clusterSizes