# Voids: 0
# Ions : 1
# Ryds : 2

function penningStartGrid(L, penningFrac)
    vect = zeros(L^2)
    rydCount = round(Int, L^2*(1-penningFrac))
    ionCount = round(Int, 0.5*L^2*penningFrac)

    vect[1:rydCount] .= 2.
    vect[rydCount+1:rydCount+ionCount] .= 1.

    grid = reshape(vect, (L,L))
    shuffle!(grid)
    return grid
end
export penningStartGrid

# Dissipates the grid in-place
function dissipateGrid!(grid, γ)
    ionPos = findall(grid .== 1.)
    switchPos = randsubseq(ionPos, γ) #Choose elements from ionPos, each with probability γ
    grid[switchPos] .= 0.
    return grid
end

#Assumes that bufferGrid is currently identical to grid
function smillaPercolationStep!(grid, bufferGrid, time, Δt, kion, k2br, γ)
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

#Based on the algorithm described by Ed that shouldn't depend on the arrangement of cells for the rate order. 
function edPercolationStep!(grid, bufferGrid, time, Δt, kion, k2br, γ)
    # Perform ionization and recombination independently
    cellCount = prod(size(grid))
    ρI = count(==(1), grid) / cellCount
    ρR = count(==(2), grid) / cellCount

    ionizationCount = round(Int, cellCount*kion*ρI*ρR*Δt)
    recombinationCount = round(Int, cellCount*k2br*ρI^2*Δt)

    # TODO: choose cells based on number of neighbours
    # Ionization
    neighbouredRs = findNeighboured(grid, 1., ==(2)) # Find all the Rydbergs that are neighboured by at least one plasma
    if length(neighbouredRs) <= ionizationCount # Less rydbergs than the number we want to ionize
        ionizationIndices = neighbouredRs #take all of them
    else
        ionizationIndices = shuffle!(neighbouredRs)[1:ionizationCount] # Selects N values without repetition
    end
    bufferGrid[ionizationIndices] .= 1.
    
    # Recombination
    neighbouredIs = findNeighboured(grid, 1., ==(1)) # Find all the plasmas that are neighboured by at least one plasma
    if length(neighbouredIs) <= recombinationCount # Less ions than the number we want to recombine
        recombinationIndices = neighbouredIs #take all of them
    else
        recombinationIndices = shuffle!(neighbouredIs)[1:recombinationCount] # Selects N values without repetition
    end
    bufferGrid[recombinationIndices] .= 2.

    grid .= bufferGrid

    # Then perform dissipation
    dissipateGrid!(grid, γ)
end

function percolate(L, dim, startGridProvider, percolationStep!, iterations, Δt, saveFract=0.1, verbose=true)
    verbose && println("Percolating the grid for $iterations iterations")
    ρRs = zeros(iterations)
    ρIs = zeros(iterations)
    size = ntuple(Returns(L), dim)
    cellCount = prod(size) # == L^dim

    saveCount = round(Int, saveFract * iterations)
    verbose && println("Saving $saveCount slices")
    saveTimes = round.(Int, range(1, iterations, saveCount))
    verbose && println("  at iterations $saveTimes")
    savedGrids = zeros(size..., saveCount)

    # Start grid, i=1
    grid = startGridProvider(L)
    @assert Base.size(grid) == size
    bufferGrid = copy(grid)
    ρRs[1] = count(==(2), grid)/cellCount
    ρIs[1] = count(==(1), grid)/cellCount
    selectdim(savedGrids, dim+1, 1) .= grid # dim+1 corresponds to the time dimension

    savei = 2
    for i in 2:iterations
        time = (i-1)*Δt
        #verbose && println("In Iteration $i")
        percolationStep!(grid, bufferGrid, time, Δt)
        copy!(bufferGrid, grid)

        ρRs[i] = count(==(2), grid)/cellCount
        ρIs[i] = count(==(1), grid)/cellCount

        if i in saveTimes
            verbose && println("Saving at time $i")
            selectdim(savedGrids, dim+1, savei) .= grid
            savei += 1
        end
    end
    return ρIs, ρRs, saveTimes, savedGrids
end
export percolate

function smillaPercolation(L, dim, penningFrac, kion, k2br, γ, iterations, saveFract = 0.1, verbose = true)
    return percolate(L, dim, 
        (l) -> penningStartGrid(l, penningFrac),
        (grid, buffer, time, Δt) -> smillaPercolationStep!(grid, buffer, time, Δt, kion, k2br, γ),
        iterations, 1., saveFract,
        verbose
    )
end
export smillaPercolation

function edPercolation(L, dim, penningFrac, kion, k2br, γ, iterations, Δt, saveFract = 0.1, verbose = true)
    return percolate(L, dim, 
        (l) -> penningStartGrid(l, penningFrac),
        (grid, buffer, time, δt) -> edPercolationStep!(grid, buffer, time, δt, kion, k2br, γ),
        iterations, Δt, saveFract,
        verbose
    )
end
export edPercolation