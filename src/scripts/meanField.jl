using DifferentialEquations
using Plots

function plasmaModel!(du, u, p, t)
    R, I, V = u
    k2br, kion, γ = p
#    dR = k2br*I^2 - kion*max(I, 0)^(max(1, 2-I))*R
    dR = k2br*I^2 - kion*I*R
    dV = γ*I
    dI = -dR-dV
    du[1] = dR
    du[2] = dI
    du[3] = dV
end
export plasmaModel!

function phaseSpace(f, u0s, p, tmax)
    sols = similar(u0s, Any)
    tspan = (0, tmax)
    for (ind, u0) in pairs(u0s)
        prob = ODEProblem(f, collect(u0), tspan, p)
        sols[ind] = solve(prob)
    end
    Rs = plot(title = "Rydberg Concentration over Time")
    Is = plot(title = "Ion Concentration over Time")
    phase = plot(title = "Phase Space", xlabel = "[NO**]", ylabel = "[NO+]")
    for sol in sols
        plot!(Rs, sol.t, sol[1,:], legend = false)
        plot!(Is, sol.t, sol[2,:], legend = false)
        plot!(phase, sol[1,:], sol[2,:], legend = false)
    end
    display(Rs)
    display(Is)
    display(phase)
    return 
end
export phaseSpace
