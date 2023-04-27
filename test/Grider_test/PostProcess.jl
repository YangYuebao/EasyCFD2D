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

function showFlow(x_uv,y_uv,uc,vc,k)
    n,m=size(x_uv)
    plt=quiver(x_uv,y_uv,quiver=(k*uc,k*vc))
    for i in [1,n]
        plot!(plt,x_uv[i,:],y_uv[i,:])
    end
    for j in [1,m]
        plot!(plt,x_uv[:,j],y_uv[:,j])
    end
    return plt
end

function showBound(bd::bound,m)
    temp = bd.fun.(range(bd.span..., length=m))
    x=Vector{Float64}(undef,0)
    y=Vector{Float64}(undef,0)
    for t in temp
        push!(x,t[1])
        push!(y,t[2])
    end
    return plot(x,y)
end



