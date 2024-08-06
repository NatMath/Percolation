using FFTW

function spatialCorrelation(grid::Matrix, cutoffDist = size(grid, 1) - 1)
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

# Much faster base time, and looks to be ~O(n^2) instead of O(n^4)
# Not sure on the correctness
# Implicitely assumes periodic boundary conditions (I think)
function spatialCorrelationFFT(grid::Matrix)
    grid .-= mean(grid)
    FT = fft(grid)
    autocorr = real.(ifft(FT .* conj.(FT))) # the i,j element contains the autocorrelation with offset i,j

    cutoff = minimum(size(grid))
    sums = zeros(cutoff^2+1)
    counts = zeros(cutoff^2+1)

    for idx in CartesianIndices(autocorr)
        d2 = sum([idx[i]^2 for i in 1:ndims(autocorr)])
        d2 <= cutoff^2+1 || continue
        sums[d2] += autocorr[idx] /prod(size(grid))
        counts[d2] += 1
    end
    corrs = sums ./ counts
    distances = sqrt.(0:cutoff^2)
    return distances[counts .!= 0], corrs[counts .!= 0]
end
export spatialCorrelationFFT

function spatialCorrelationFFT(grid)
    avg = mean(grid)
    FT = fft(grid .- avg) #Real fft, optimized compared to the most generic fft
    convolution = real.(ifft(FT .* conj.(FT)))

    maxd2 = sum(size(grid) .^2 ) # maximum distance is the opposite corner
    sums = zeros(maxd2+1)
    counts = zeros(maxd2+1)

    for idx in CartesianIndices(grid)
        d2 = sum([idx[i]^2 for i in 1:ndims(grid)])
        sums[d2] += convolution[idx]
        counts[d2] += 1
    end
    corrs = sums ./ (counts .* prod(size(grid)))
    distances = sqrt.(0:maxd2)
    return distances[counts .!= 0], corrs[counts .!= 0]
end