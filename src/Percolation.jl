module Percolation
using Revise
using Random
using Plots
using Dates

include("scripts/gridDensity.jl")
include("scripts/meanField.jl")
include("scripts/SFIDistribution.jl")
include("scripts/statsTest.jl")

include("tools/boundaryConditions.jl")
include("tools/gridTools.jl")
include("tools/plotTools.jl")
include("tools/mathTools.jl")
include("tools/statsTools.jl")

include("steppedPercolation.jl")

include("analytics/singlePercolation.jl")
include("analytics/rateOrders.jl")
include("analytics/edSmillaCompare.jl")


end # module Percolation
