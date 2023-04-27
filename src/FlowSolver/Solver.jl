function solvefield(coodinate::AbstractCoodinateTypes, convectionSecheme::T, x_uv::Matrix{Float64}, y_uv::Matrix{Float64}, mu::Float64, rho::Float64,bounds::Vector{bound}; abstol::Float64=1e-5, maxiter::Int64=100) where T<:AbstractScheme
    n,m=size(x_uv)
    u, v, p, x_u, y_u, x_v, y_v, alpha, beta, gamma, Ja = fieldA(x_uv, y_uv)
    row, col = index_generation(n, m, getPlace(convectionSecheme))
    rowp, colp = index_generation(n, m, getPlace(Point5()))
    err::Float64 = Inf
    count::Int64 = 0
    #=
        mu = 1.0
        rho = 1.0
    =#
    while err > abstol && count < maxiter
        valu = Vector{Float64}(undef, 0)
        bu = Vector{Float64}(undef, 0)
        valv = Vector{Float64}(undef, 0)
        bv = Vector{Float64}(undef, 0)
        Apu = zeros(n, m)
        Apv = zeros(n, m)
        U, V = getUV(u, v, x_u, x_v, y_u, y_v)
        p_u, p_v = fieldDiff(p)

        Su = p_v .* y_u - p_u .* y_v
        Sv = p_u .* x_v - p_v .* x_u
        # 生成系数
        renew_coff_field!(coodinate,convectionSecheme, n, m, mu, rho, valu, bu, alpha, beta, gamma, Ja,x_u,x_v,y_u,y_v, U, V, u, v, x_uv, y_uv, u, Su, Apu, :u, bounds)

        renew_coff_field!(coodinate,convectionSecheme, n, m, mu, rho, valv, bv, alpha, beta, gamma, Ja,x_u,x_v,y_u,y_v, U, V, u, v, x_uv, y_uv, v, Sv, Apv, :v, bounds)

        A = sparse(row, col, valu)
        u = reshape(A \ bu, n, m)
        A = sparse(row, col, valv)
        v = reshape(A \ bv, n, m)


        valp = Vector{Float64}(undef, 0)
        bp = Vector{Float64}(undef, 0)

        SIMPLE(coodinate, n, m, valp, bp, rho, x_u,x_v,y_u,y_v, Apu, Apv, p, U, V,u,v,x_uv,y_uv,Ja, alpha, gamma, :p, bounds)

        A = sparse(rowp, colp, valp)

        p += reshape(A \ bp, n, m)
        
        count+=1
        err=maximum(abs.(bp))
    end
    return u,v,p
end