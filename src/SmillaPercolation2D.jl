# Voids: 0
# Ions : 1
# Ryds : 2

function startGrid(L, penningFrac)
    vect = zeros(L^2)
    rydCount = round(Int, L^2*(1-penningFrac))
    ionCount = round(Int, 0.5*L^2*penningFrac)

    vect[1:rydCount] .= 2.
    vect[rydCount+1:rydCount+ionCount] .= 1.

    grid = reshape(vect, (L,L))
    shuffle!(grid)
    return grid
end
export startGrid

# Dissipates the grid in-place
function dissipateGrid!(grid, γ)
    ionPos = findall(grid .== 1.)
    switchPos = randsubseq(ionPos, γ) #Choose elements from ionPos, each with probability γ
    grid[switchPos] .= 0.
    return grid
end
export dissipateGrid!

#Assumes that bufferGrid is currently identical to grid
function percolationStep!(grid, kion, k2br, γ, bufferGrid = copy(grid))
    # Perform ionization and recombination independently
    ρI = count(==(1), grid) / prod(size(grid))
    ρR = count(==(2), grid) / prod(size(grid))

    # Ionization
    neighbouredRs = findNeighboured(grid, 1., ==(2)) # Find all the Rydbergs that are neighboured by at least one plasma
    ionizationp = clamp(kion * ρI, 0., 1.)
    ionizationIndices = randsubseq(neighbouredRs, ionizationp)
    bufferGrid[ionizationIndices] .= 1.
    
    # Recombination
    neighbouredIs = findNeighboured(grid, 1., ==(1)) # Find all the plasmas that are neighboured by at least one plasma
    recombinationp = clamp(k2br * ρI, 0., 1.)
    recombinationIndices = randsubseq(neighbouredIs, recombinationp)
    bufferGrid[recombinationIndices] .= 2.

    grid .= bufferGrid

    # Then perform dissipation
    dissipateGrid!(grid, γ)
end
export percolationStep!

function percolate(L, penningFrac, kion, k2br, γ, iterations, saveFract = 0.1, verbose = true)
    verbose && println("Percolating the grid for $iterations iterations")
    ρRs = zeros(iterations)
    ρIs = zeros(iterations)
    cellCount = L^2

    saveCount = round(Int, saveFract * iterations)
    verbose && println("Saving $saveCount slices")
    saveTimes = round.(Int, range(1, iterations, saveCount))
    verbose && println("  at iterations $saveTimes")
    savedGrids = zeros(L, L, saveCount)

    # Start grid, i=1
    grid = startGrid(L, penningFrac)
    bufferGrid = copy(grid)
    ρRs[1] = count(==(2), grid)/cellCount
    ρIs[1] = count(==(1), grid)/cellCount
    savedGrids[:,:,1] .= grid

    savei = 2
    for i in 2:iterations
        #verbose && println("In Iteration $i")
        percolationStep!(grid, kion, k2br, γ, bufferGrid)
        copy!(bufferGrid, grid)

        ρRs[i] = count(==(2), grid)/cellCount
        ρIs[i] = count(==(1), grid)/cellCount

        if i in saveTimes
            verbose && println("Saving at time $i")
            savedGrids[:,:,savei] .= grid
            savei += 1
        end
    end
    return ρIs, ρRs, saveTimes, savedGrids
end
export percolate