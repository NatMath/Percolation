

function singlePercolation()
    # Simulation Parameters
    L = 200
    penningFrac = 0.2
    kion = k2br = 1.
    γ = kion/1e2

    itr = Int(1e3)
    
    ρIs, ρRs, saveTimes, savedGrids = percolate(L, penningFrac, kion, k2br, γ, itr, 0.01, true)

    saveCounts = length(saveTimes)
 

    #Plotting
    plotlyjs()

    ρIsPlot = plot(ρIs, title = "Ion density", xlabel = "Time (itr)", ylabel = "Density (/site)")
    display(ρIsPlot)

    ρRsPlot = plot(ρRs, title = "Rydberg density", xlabel = "Time (itr)", ylabel = "Density (/site)")
    display(ρRsPlot)

    phasePlot = plot(ρRs, ρIs, title = "Phase Space", xlabel = "Rydberg density (/site)", ylabel = "Ion density (/site)")
    display(phasePlot)

    # Autocorrelations
    autoCorrIs = round.(Int, range(1, saveCounts, 10))
    print("Calculating autocorrelations at indices $autoCorrIs")

    gr()
    autoCorrPlot = plot(1:10, 1:10)
    autoCorrAnim = @animate for i in autoCorrIs
        grid = savedGrids[:,:,i]
        t = saveTimes[i]
        distances, autoCorrs = spatialCorrelation(grid .== 1.) # Autocorrelation of ions
        distances, autoCorrs = recMeanSmooth(distances, autoCorrs, 2)
        autoCorrs = max.(abs.(autoCorrs), 1e-5)
        plot(distances[2:end], autoCorrs[2:end], xaxis = :log, yaxis = :log, title = "Time $t")
    end
    name = timestamp() * " - " * "Autocorrelation Plot"
    path = "./plots"
    println("Saving file at: ", path)
    println("under name: ", name)
    fullpath = path * "/" * name * ".gif"
    gif(autoCorrAnim, fullpath, fps=1)


end
export singlePercolation