function timestamp()
    return Dates.format(now(), "YYYY_mm_dd__HH_MM_SS")
end
export timestamp


function plotGrid(grid)
    plotlyjs()
    L = size(grid)[1]
    colors = palette([:black, :red, :blue])
    heatmap(grid, c=colors, legend = :none, xlims=(1,L), ylims = (1,L), aspect_ratio=:equal)
end
export plotGrid

function animateGrids(grids::Array{<:Number, 3}, timestamps = 1:size(grids,3);
         path = "./plots", name = "Grid Animation", fps=15)
    gr()
    L = size(grids,1)
    colors = palette([:black, :red, :blue])
    hmap = heatmap(grids[:,:,begin], title = timestamps[1], c=colors, legend = :none, xlims=(1,L), ylims = (1,L), aspect_ratio=:equal)
    anim = @animate for i ∈ 1:size(grids, 3)
        heatmap!(hmap[1], grids[:,:,i], title = timestamps[i], c=colors, legend = :none, xlims=(1,L), ylims = (1,L), aspect_ratio=:equal)
    end
    name = timestamp() * " - " * name
    println("Saving file at: ", path)
    println("under name: ", name)
    fullpath = path * "/" * name * ".gif"
    gif(anim, fullpath, fps=fps)
end
export animateGrids