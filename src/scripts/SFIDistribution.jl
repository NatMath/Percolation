using Plots 

function plaw(α, xmin = 1., xmax = 1e6) #Returns a random sample of the power-law between xmin and infinity
    r = rand() #Uniformly
    sample = xmin * (1-r) ^ (-1/(α-1))
    return sample
end
export plaw

function binData(data, nbins = 30, xmin = 1.)
    xmax = maximum(data)
    binStarts = 10 .^ range(log10(xmin), log10(xmax), nbins+1)
    binCenters = (binStarts[begin:end-1] + binStarts[begin+1:end])/2

    counts = [
        count(x -> binStarts[i] <= x < binStarts[i+1], data) for i in 1:nbins
    ]
    return binCenters, counts
end
export binData

function testDist(n, α=1.4, nbins = n/100)
    data = [plaw(α, 1.) for i in 1:n]
    bins, counts = binData(data, nbins, 1.)

    nonZeroBins = bins[counts .!= 0]
    nonZeroCounts = counts[counts .!= 0]
    distPlot = plot(nonZeroBins, nonZeroCounts/length(data),
         xaxis=:log, yaxis=:log, title = "P-law Probability Distribution", xlabel = "N", ylabel = "p(N)")
    display(distPlot)

    cummCounts = [sum(counts[i:end]) for i in 1:nbins]
    cummDistPlot = plot(bins[cummCounts .!= 0], cummCounts[cummCounts .!= 0]/length(data),
        xaxis=:log, yaxis=:log, title = "P-law Probability Cummulative Distribution", xlabel = "N", ylabel = "P'(N)")
    display(cummDistPlot)

    nSamples = Int(n/1e3)
    nSums = Int(1e7)
    sampledSums = [sum(rand(data, nSamples)) for i in 1:nSums]
    #println(sampledSums)
    sumBins, sumCounts = binData(sampledSums, 100, minimum(sampledSums))
    plot(sumBins[sumCounts .!= 0], sumCounts[sumCounts .!= 0],
         xaxis=:log, yaxis = :log, xlabel = "Sum samples", title = "Sample sum histogram")

end
export testDist
