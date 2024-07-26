# Sample from plaw distribution, test binning/plotting
function plawSampling()
    N = 100000 #Number of samples to take
    α = 1.4
    shape = α-1
    dist = Pareto(shape)
    samples = rand(dist, N)
    
    Nbins = round(Int, N/100)
    bins, counts = trimZeros(logBinning(samples, Nbins, min = 1)...)
    scatter(bins, counts, xaxis=:log, yaxis=:log)
end
export plawSampling

function sumOfSamples()
    m = 10000
    N = 10000
    α=1.4
    shape=α-1
    dist = Pareto(shape)
    samples = [sum(rand(dist, N)) for _ in 1:m]

    Nbins = round(Int, N/100)
    bins, counts = trimZeros(logBinning(samples, Nbins, min = 1)...)
    scatter(bins, counts, xaxis=:log, yaxis=:log, label="Sum of plaw")

    μ, σ = 2e7, 1e9
    scale = 1e6
    levySamples = rand(Levy(μ, σ), 100000)
    plot!(trimZeros(logBinning(levySamples, 1000, min=1)...)...,xaxis=:log, yaxis=:log, label = "Sampled Lévy")
    #xs = LinRange(minimum(samples), maximum(samples), 1000)
    #ys = @. scale*sqrt(σ/(2π*(xs-μ)^3))*exp(-σ/(2*(xs-μ)))
    #plot!(xs, scale * ys, label = "Lévy PDF")
end
export sumOfSamples

function plotPDFs()
    xs = 10 .^ LinRange(0, 4, 100)

    α=1.4
    ys_plaw = @. (α-1)*xs^(-α)

    μ, σ = 2, .25
    ys_lognormal = @. 1/(xs*σ*sqrt(2π))*exp(-((log(xs)-μ)^2)/(2σ^2))

    μ, σ = 1, 1.
    ys_levy = @. sqrt(σ/(2π*(xs-μ)^3))*exp(-σ/(2*(xs-μ)))

    plot(trimZeros(xs, ys_plaw)..., xaxis=:log, yaxis=:log, label = "plaw")
    #plot!(trimZeros(xs, ys_lognormal)..., label = "Lognormal")
    plot!(trimZeros(xs, ys_levy)..., label = "Levy")

    #Log-Normal is significantly sub-plaw. Heavy-tailed, but not by much
end
export plotPDFs
