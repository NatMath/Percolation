function edSmillaCompare()
    # Simulation Parameters
    L = 200
    penningFrac = 0.2
    kion = k2br = 1.
    γ = kion/1e2
    Δt = .1

    itr = Int(1e4)
    
    ρIs_S, ρRs_S, saveTimes_S, savedGrids_S = smillaPercolation(L, 2, penningFrac, kion*Δt, k2br*Δt, γ*Δt, itr, 0.01, true)
    ρIs_E, ρRs_E, saveTimes_E, savedGrids_E = edPercolation(    L, 2, penningFrac, kion,    k2br,    γ*Δt, itr, Δt, 0.01, true)

    saveCounts = length(saveTimes_S)
 

    #Plotting
    plotlyjs()

    ρIsPlot = plot(ρIs_S, title = "Ion density", xlabel = "Time (itr)", ylabel = "Density (/site)", label="Smilla")
    plot!(ρIsPlot, ρIs_E, title = "Ion density", xlabel = "Time (itr)", ylabel = "Density (/site)", label="Ed")
    display(ρIsPlot)

    ρRsPlot = plot(ρRs_S, title = "Rydberg density", xlabel = "Time (itr)", ylabel = "Density (/site)", label="Smilla")
    plot!(ρRsPlot, ρRs_E, title = "Rydberg density", xlabel = "Time (itr)", ylabel = "Density (/site)", label="Ed")
    display(ρRsPlot)

    phasePlot = plot(ρRs_S, ρIs_S, title = "Phase Space", xlabel = "Rydberg density (/site)", ylabel = "Ion density (/site)", label="Smilla")
    plot!(phasePlot, ρRs_E, ρIs_E, title = "Phase Space", xlabel = "Rydberg density (/site)", ylabel = "Ion density (/site)", label="Ed")
    display(phasePlot)

end
export edSmillaCompare