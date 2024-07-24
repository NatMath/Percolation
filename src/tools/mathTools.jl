# Returns the mean of the array along the given dimensions
mean(array) = sum(array) / prod(size(array))
function mean(array, dims)
    return sum(array, dims=dims) ./ prod(size(array)[[dims...]])
end
export mean

# Smooths a vector by taking 
function meanSmooth(xs, ys, meanDst) 
    L = length(xs)
    newxs = xs[begin+meanDst:end-meanDst]
    newys = [sum(ys[i - meanDst:i + meanDst])/(2*meanDst+1) for i in (1+meanDst):(L-meanDst)]
    return (newxs, newys)
end
export meanSmooth

function recMeanSmooth(xs, ys, itrs::Int)
    itrs == 0 && return (xs, ys)
    L = length(xs)
    println(itrs)
    println(length(xs))
    println(length(ys))
    excess = mod(L, 3)
    newxs = [xs[2:3:L-1]; xs[end-excess + 1:end]]
    newys = [[mean(ys[i-1:i+1]) for i in 2:3:L-1]; ys[end-excess + 1:end]]
    return recMeanSmooth(newxs, newys, itrs - 1)
end
export recMeanSmooth