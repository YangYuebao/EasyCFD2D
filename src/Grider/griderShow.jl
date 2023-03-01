using Plots

begin
    include("bounds3.jl")
    include("GS_solver.jl")
    plt = plot(
        legend=:none,
        xlim=(0, 3),
        ylim=(0, 2)
    )
    for i = 1:n
        plot!(plt, x_uv[i, :][:], y_uv[i, :][:], linecolor=:blue)
    end

    for j = 1:m
        plot!(plt, x_uv[:, j][:], y_uv[:, j][:], linecolor=:blue)
    end
    savefig(plt, "./bodyFitGrid/plt1.png")
end


begin
    include("bounds2.jl")
    include("GS_solver.jl")
    plt = plot(
        legend=:none,
        xlim=(-1.5, 1.5),
        ylim=(0, 2)
    )

    for i = 1:m
        plot!(plt, x_uv[i, :][:], y_uv[i, :][:], linecolor=:blue)
    end

    for j = 1:n
        plot!(plt, x_uv[:, j][:], y_uv[:, j][:], linecolor=:blue)
    end
    savefig(plt, "./bodyFitGrid/plt2.png")
end
