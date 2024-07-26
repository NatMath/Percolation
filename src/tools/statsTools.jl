using Distributions

function logBinning(data, Nbins; min = minimum(data), max = maximum(data))
    logmin = log(min)
    logmax = log(max)
    logedges = LinRange(logmin, logmax, Nbins+1)
    edges = exp.(logedges)
    logcentres = (logedges[begin:end-1] + logedges[begin+1:end])/2
    centres = exp.(logcentres)

    counts = [
        count(@. (start <= data) & (data < stop)) for (start, stop) in zip(edges[begin:end-1], edges[begin+1:end])
    ]
    return centres, counts
end
export logBinning

function linBinning(data, Nbins; min = minimum(data), max = maximum(data))
    edges = LinRange(min, max, Nbins+1)
    centres = (edges[begin:end-1] + edges[begin+1:end])/2

    counts = [
        count(@. start <= data & data < stop) for (start, stop) in zip(edges[begin:end-1], edges[begin+1:end])
    ]
    return centres, counts
end
export linBinning

function trimZeros(xs, ys) # Returns all xs and ys where y =/= 0
    newXs = xs[ys .!= 0]
    newYs = ys[ys .!= 0]
    return newXs, newYs
end
