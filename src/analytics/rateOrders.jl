using CurveFit
using Interpolations

function rateOrderSim(autoAnalyse=true)
    mode = (:ion,) # :ion, :recomb
    ρs = 0:0.01:0.75
    nρs = length(ρs)

    samples = 400

    L = 100

    kion = 1.
    k2br = 1.

    rateIon = zeros(nρs, nρs)
    if :ion in mode #Ionization
        println("Ionization")
        vect = zeros(L^2)
        bufferGrid = zeros((L,L))
        for i in 1:nρs, j in 1:nρs
            ρI = ρs[i]
            ρR = ρs[j]
            rydCount = round(Int, ρR * L^2)
            ionCount = round(Int, ρI * L^2)
            if rydCount + ionCount > L^2
                continue
            end
            for iter in 1:samples
                vect .= 0
                vect[1:rydCount] .= 2.0
                vect[rydCount+1:rydCount+ionCount] .= 1.0
            
                grid = reshape(vect, (L, L))
                shuffle!(grid)
                copy!(bufferGrid, grid)
            
                percolationStep!(grid, kion, 0.0, 0.0, bufferGrid)
                rateIon[j, i] += -(count(==(2), grid) - rydCount) / L^2
            end
            if j == 1
                println("Outer loop iteration $i/$nρs done")
            end
        end
        rateIon ./= samples
    end

    rateRecomb = zeros(nρs, nρs)
    if :recomb in mode #2br
        println("Recombination")
        vect = zeros(L^2)
        bufferGrid = zeros((L,L))
        for i in 1:nρs, j in 1:nρs
            ρI = ρs[i]
            ρR = ρs[j]
            rydCount = round(Int, ρR * L^2)
            ionCount = round(Int, ρI * L^2)
            if rydCount + ionCount > L^2
                continue
            end
            for iter in 1:samples
                vect .= 0
                vect[1:rydCount] .= 2.0
                vect[rydCount+1:rydCount+ionCount] .= 1.0
            
                grid = reshape(vect, (L, L))
                shuffle!(grid)
                copy!(bufferGrid, grid)
            
                percolationStep!(grid, 0.0, k2br, 0.0, bufferGrid)
                rateRecomb[j, i] += (count(==(2), grid) - rydCount) / L^2
            end
            if j == 1
                println("Outer loop iteration $i/$nρs done")
            end
        end
        rateRecomb ./= samples
    end  
    autoAnalyse && rateOrderAnalysis(mode, ρs, ρs, rateIon, rateRecomb, kion, k2br)
    return (mode, ρs, ρs, rateIon, rateRecomb, kion, k2br)
end
export rateOrderSim

function rateOrderAnalysis(mode, ρRs, ρIs, rateIon, rateRecomb, kion, k2br)
    nρRs = length(ρRs)
    nρIs = length(ρIs)

    plotlyjs()

    #===========================================================================================================
                            Ionizations
    ===========================================================================================================#
    if :ion in mode
        rydSliceIndex = 10#round(Int, nρs/10)
        rydSlice = rateIon[:,rydSliceIndex]
        rydOrder = power_fit(ρRs[rydSlice .!= 0], rydSlice[rydSlice .!= 0])
        println("P-fit params for Ionization rydberg dependence (mul, power): $rydOrder")
        rydInterp = linear_interpolation(ρRs, rydSlice)
        initialRydOrder = power_fit((ρRs[rydSlice .!= 0])[begin:round(Int, nρRs/10)],
                                    (rydSlice[rydSlice .!= 0])[begin:round(Int, nρRs/10)])[2]
        #diffeqRateOrder(rydInterp, (ρRs[rydSlice .!= 0])[[1,end]], initialRydOrder,
        #     toplot=true, title = "Rate order in Ryd for ionization", xlabel = "Rydberg Density", ylabel = "Rate Order")
        rydOrders = [rateOrderLaurent(rydInterp, ρ, 4, iters = 20, atol=1e-6, verbose=false) for ρ=ρRs]
        display(plot(meanSmooth(ρRs, rydOrders, 3)..., xlabel = "Rydberg density", title = "Ionization order in Rydberg"))

        ionSliceIndex = 10
        ionSlice = rateIon[ionSliceIndex,:]
        ionOrder = power_fit(ρIs[ionSlice .!= 0], ionSlice[ionSlice .!= 0])
        println("P-fit params for Ionization ion dependence (mul, power): $ionOrder")
        ionInterp = linear_interpolation(ρIs, ionSlice)
        initialIonOrder = power_fit((ρIs[ionSlice .!= 0])[begin:round(Int, nρIs/10)],
                                    (ionSlice[ionSlice .!= 0])[begin:round(Int, nρRs/10)])[2]
        #diffeqRateOrder(ionInterp, (ρIs[ionSlice .!= 0])[[1,end]], initialIonOrder,
        #     toplot=true, title = "Rate order in Ion for ionization", xlabel = "Ion Density", ylabel = "Rate Order")
        ionOrders = [rateOrderLaurent(ionInterp, ρ, 4, iters = 20, atol=1e-6, verbose=false) for ρ=ρIs]
        display(plot(meanSmooth(ρIs, ionOrders, 3)..., xlabel = "Ion density", title = "Ionization order in Ion"))

        pRyd = plot(ρRs, rydSlice, 
            xlabel = "Rydberg Density", title = "Ionization Reaction Rate, kion = $kion \n Ion density = $(ρIs[rydSliceIndex])")
        plot!(pRyd, ρRs, rydOrder[1] .* ρRs .^ rydOrder[2])
        display(pRyd)
        pIon = plot(ρIs, ionSlice, 
            xlabel = "Ion Density", title = "Ionization Reaction Rate, kion = $kion \n Rydberg density = $(ρRs[ionSliceIndex])")
        plot!(pIon, ρIs, ionOrder[1] .* ρIs .^ ionOrder[2])
        display(pIon)
        display(plot(ρRs, ρIs, rateIon, st=:surface, camera = (-30, 30),
            xlabel = "Rydberg Density", ylabel = "Ion Density", title = "Ionization Reaction Rate, kion=$kion"))
    end

    #===========================================================================================================
                            Recombination
    ===========================================================================================================#
    if :recomb in mode
        rydSliceIndex = 10#round(Int, nρs/10)
        rydSlice = rateRecomb[:,rydSliceIndex]
        rydOrder = power_fit(ρRs[rydSlice .!= 0], rydSlice[rydSlice .!= 0])
        println("P-fit params for Recombination rydberg dependence (mul, power): $rydOrder")
        rydInterp = linear_interpolation(ρRs, rydSlice)
        initialRydOrder = power_fit((ρRs[rydSlice .!= 0])[begin:round(Int, nρRs/10)],
                                    (rydSlice[rydSlice .!= 0])[begin:round(Int, nρRs/10)])[2]
        #diffeqRateOrder(rydInterp, (ρRs[rydSlice .!= 0])[[1,end]], initialRydOrder,
        #     toplot=true, title = "Rate order in Ryd for ionization", xlabel = "Rydberg Density", ylabel = "Rate Order")
        rydOrders = [rateOrderLaurent(rydInterp, ρ, 4, iters = 20, atol=1e-6, verbose=false) for ρ=ρRs]
        display(plot(ρRs, rydOrders, xlabel = "Rydberg density", title = "Recombination order in Rydberg"))

        ionSliceIndex = 10
        ionSlice = rateRecomb[ionSliceIndex,:]
        ionOrder = power_fit(ρIs[ionSlice .!= 0], ionSlice[ionSlice .!= 0])
        println("P-fit params for Recombination ion dependence (mul, power): $ionOrder")
        ionInterp = linear_interpolation(ρIs, ionSlice)
        initialIonOrder = power_fit((ρIs[ionSlice .!= 0])[begin:round(Int, nρIs/10)],
                                    (ionSlice[ionSlice .!= 0])[begin:round(Int, nρRs/10)])[2]
        #diffeqRateOrder(ionInterp, (ρIs[ionSlice .!= 0])[[1,end]], initialIonOrder,
        #     toplot=true, title = "Rate order in Ion for ionization", xlabel = "Ion Density", ylabel = "Rate Order")
        ionOrders = [rateOrderLaurent(ionInterp, ρ, 4, iters = 20, atol=1e-6, verbose=false) for ρ=ρIs]
        display(plot(ρIs, ionOrders, xlabel = "Ion density", title = "Recombination order in Ion"))

        pRyd = plot(ρRs, rydSlice, 
            xlabel = "Rydberg Density", title = "Recombination Reaction Rate, k2br = $k2br \n Ion density = $(ρIs[rydSliceIndex])")
        plot!(pRyd, ρRs, rydOrder[1] .* ρRs .^ rydOrder[2])
        display(pRyd)
        pIon = plot(ρIs, ionSlice, 
            xlabel = "Ion Density", title = "Recombination Reaction Rate, k2br = $k2br \n Rydberg density = $(ρRs[ionSliceIndex])")
        plot!(pIon, ρIs, ionOrder[1] .* ρIs .^ ionOrder[2])
        display(pIon)
        display(plot(ρRs, ρIs, rateRecomb, st=:surface, camera = (-30, 30),
            xlabel = "Rydberg Density", ylabel = "Ion Density", title = "Recombination Reaction Rate, k2br=$k2br"))
    end
    return (mode, ρRs, ρIs, rateIon, rateRecomb, kion, k2br)
end
export rateOrderAnalysis

function diffeq(n, args, x) # Assuming f=x^n, knowing f and wanting to find n
    f = args
    return (x/f(x)*ForwardDiff.derivative(f, x) - n)/(x*log(x))
end
function diffeqRateOrder(f, xs, g0; verbose=false, toplot=false, kwargs...)
    verbose && println("Initial value of diffeq: $(diffeq(g0, Nothing, xs[1]))")
    prob = ODEProblem(diffeq, g0, xs, f)
    sol = solve(prob, reltol = 1e-6)
    if toplot
        p = plot(sol.t, sol[:]; kwargs...)
        display(p)
    end
end
export diffeqRateOrder

























# Rate order of ys = f.(xs)
# Assumes that the order is constant, otherwise doesn't work
# Discrete difference method for derivative
function rateOrderFit(xs, ys)
    #return xs[2:end] .* diff(log.(ys)) ./ diff(log.(xs))
    return power_fit(xs, ys)
end
export rateOrder


# Me being dumb for 3h. Some cool diffeqs in there, but can be solved with just a log...
using ForwardDiff
using DifferentialEquations

function rateOrderLaurent(f, x::Real, _maxn; iters = 4, atol = 1e-1, verbose=true)
#=     println("A")
    derivs = Vector{Any}(undef, nmax)
    derivs[1] = x -> ForwardDiff.derivative(f, x)
    for i in 2:nmax
        derivs[i] = x -> ForwardDiff.derivative(derivs[i-1], x)
    end

    derivsAtX = [derivs[i](x) for i in 1:nmax]
    println(derivsAtX)

    maxn = findfirst(derivsAtX .<= 0) - 1
    minn = maxn - 1
    println("Order of function is between $minn and $maxn")
    println("Derivs of f/x^minn: $(ForwardDiff.derivative(x -> f(x)/x^minn, x))")
    println("Derivs of f/x^maxn: $(ForwardDiff.derivative(x -> f(x)/x^maxn, x))")

    if derivsAtX[maxn+1] == 0.
        println("Exact integer order $maxn")
        return maxn
    end
 =#

 
    minn = 0.
    maxn = convert(Float64, _maxn)
    for i in 1:iters
        testn = (minn + maxn)/2
        verbose && println("Testing at n=$testn")
        deriv = ForwardDiff.derivative(x -> f(x)/x^testn, x)
        verbose && println("Deriv = $deriv")
        if abs(deriv) <= atol/x^testn
            verbose && println("Found n, terminating at iteration $i")
            return testn
        elseif deriv > 0
            minn = testn
        else
            maxn = testn
        end
    end
    verbose && println("Didn't converge in $iters iterations with tol $atol")
    verbose && println("Order of function is between $minn and $maxn")
    verbose && println("Derivs of f/x^minn: $(ForwardDiff.derivative(x -> f(x)/x^minn, x))")
    verbose && println("Derivs of f/x^maxn: $(ForwardDiff.derivative(x -> f(x)/x^maxn, x))")
    return (minn+maxn)/2
end
export rateOrderLaurent

function rateOrder(f, xs::Vector{<:Real}, _maxn; iters=4, atol=1e-3)
    #orders = [rateOrder(f, x, _maxn, iters=iters, atol=atol, verbose=false) for x in xs]
    orders = @. log(f(xs))/log(xs)
    p = plot(xs, orders)
    #plot!(p, xs, @. 1+xs/100 + xs*log(xs)*0.01)
    display(p)
    return orders
end

