using Plots
function gridPlot(x_uv, y_uv)
    n,m=size(x_uv)
    plt = plot(
        legend=:none,
        aspect_ratio=:equal
    )
    for i = 1:n
        plot!(plt, x_uv[i, :][:], y_uv[i, :][:], linecolor=:blue)
    end

    for j = 1:m
        plot!(plt, x_uv[:, j][:], y_uv[:, j][:], linecolor=:blue)
    end
    return plt
end