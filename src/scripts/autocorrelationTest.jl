
function autocorrelationTest()
    L = 200

    # 2D testing
    noise = rand(L, L)
    per1d = cos.(collect(1:L) ./ 5) .* ones(200)'
    per2d = cos.(collect(1:L) ./ 5) .* cos.(collect(1:L) ./ 5)'
    #display(heatmap(per2d))

    noisePlot = plot(spatialCorrelation(noise)..., label = "Original")
    plot!(noisePlot, spatialCorrelationFFT(noise)..., label = "FFT")
    display(noisePlot)

    per1dPlot = plot(spatialCorrelation(per1d)..., label = "Original")
    plot!(per1dPlot, spatialCorrelationFFT(per1d)..., label = "FFT")
    display(per1dPlot)

    per2dPlot = plot(spatialCorrelation(per2d)..., label = "Original")
    plot!(per2dPlot, spatialCorrelationFFT(per2d)..., label = "FFT")
    display(per2dPlot)

end
export autocorrelationTest