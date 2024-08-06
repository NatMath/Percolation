

function singlePercolation()
    # Simulation Parameters
    L = 200
    dims = 2

    penningFrac = 0.2
    kion = k2br = 1.
    γ = kion/1e2

    itr = Int(1e3)
    
    ρIs, ρRs, saveTimes, savedGrids = smillaPercolation(L, dims, penningFrac, kion, k2br, γ, itr, 0.01, true)

    saveCounts = length(saveTimes)
 

    #Plotting
    plotlyjs()

    ρIsPlot = plot(ρIs, title = "Ion density", xlabel = "Time (itr)", ylabel = "Density (/site)")
    display(ρIsPlot)

    ρRsPlot = plot(ρRs, title = "Rydberg density", xlabel = "Time (itr)", ylabel = "Density (/site)")
    display(ρRsPlot)

    phasePlot = plot(ρRs, ρIs, title = "Phase Space", xlabel = "Rydberg density (/site)", ylabel = "Ion density (/site)")
    display(phasePlot)

    # Cluster sizes
    gr()
    clusterSizePlot = plot()
    clusterSizeAnim = @animate for i in 1:saveCounts
        grid = Array(selectdim(savedGrids, dims+1, i))
        t = saveTimes[i]
        sizes = clusterSizes(grid .== 2, NormalBoundary())
        bins, counts = trimZeros(logBinning(sizes, ceil.(Int, length(sizes)/100))...)
        plot(bins, counts, xaxis=:log, yaxis=:log, title = "Time $t")
    end
    name = timestamp() * " - " * "Cluster Distribution Plot"
    path = "./plots"
    println("Saving file at: ", path)
    println("under name: ", name)
    fullpath = path * "/" * name * ".gif"
    gif(clusterSizeAnim, fullpath, fps=1)


    #=
    # Autocorrelations
    autoCorrIs = round.(Int, range(1, saveCounts))
    print("Calculating autocorrelations at indices $autoCorrIs")

    gr()
    autoCorrPlot = plot
    autoCorrAnim = @animate for i in autoCorrIs
        grid = savedGrids[:,:,i]
        t = saveTimes[i]
        distances, autoCorrs = spatialCorrelationFFT(grid .== 1.) # Autocorrelation of ions
        #distances, autoCorrs = recMeanSmooth(distances, autoCorrs, 2)
        #autoCorrs = max.(abs.(autoCorrs), 1e-5)
        autoCorrs .= abs.(autoCorrs)
        plot(distances, autoCorrs, xaxis = :log, yaxis = :log, title = "Time $t")
    end
    name = timestamp() * " - " * "Autocorrelation Plot"
    path = "./plots"
    println("Saving file at: ", path)
    println("under name: ", name)
    fullpath = path * "/" * name * ".gif"
    gif(autoCorrAnim, fullpath, fps=1)
    =#

end
export singlePercolation