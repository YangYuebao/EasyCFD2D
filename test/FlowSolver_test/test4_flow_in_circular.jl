using EasyCFD2D,SparseArrays

begin
    r1 = 1
    r2 = 2

    function bd1(t)
        [-r2 + t, 0]
    end
    t1 = [0, r2 - r1]

    #bd2 t in [-1,1]
    function bd2(t)
        r1 * [cos(t), sin(t)]
    end
    t2 = [pi, 0]

    #bd3 t in [0,pi/2]
    function bd3(t)
        [r1 + t, 0]
    end
    t3 = [0, r2 - r1]

    #bd4 t in [2,-2]
    function bd4(t)
        r2 * [cos(t), sin(t)]
    end
    t4 = [0, pi]

    bounds = Vector{bound}(undef, 4)
    bounds[1] = bound(bd1, t1, pressureInlet(), 10)
    bounds[2] = bound(bd2, t2, stillWall())
    bounds[3] = bound(bd3, t3, FDOutlet(), 0)
    bounds[4] = bound(bd4, t4, stillWall())

    include("../Grider_test/PostProcess.jl")
    #x_uv,y_uv=GSGrider(m,n,bounds)
    m = 30
    n = 30
    #x_uv,y_uv=GSGrider(m,n,bounds)
    x_uv, y_uv = EasyCFD2D.JacobianGrider(n, m, bounds, maxep=1e-5, relax=0.2, displayStep=10,maxcount=400)
    #gridPlot(x_uv,y_uv)
    convectionSecheme=SecondOrderUpwind()
    coodinate=Rectangular()
    u, v, p, x_u, y_u, x_v, y_v, alpha, beta, gamma, Ja = fieldA(x_uv, y_uv)
    row, col = EasyCFD2D.index_generation(n, m, getPlace(convectionSecheme))
    rowp, colp = EasyCFD2D.index_generation(n, m, getPlace(EasyCFD2D.Point5()))
    mu = 1.0
    rho = 1.0
    U, V = getUV(u, v, x_u, x_v, y_u, y_v)
end

begin
    valu = Vector{Float64}(undef, 0)
    bu = Vector{Float64}(undef, 0)
    valv = Vector{Float64}(undef, 0)
    bv = Vector{Float64}(undef, 0)
    Apu = zeros(n, m)
    Apv = zeros(n, m)
    p_u, p_v = fieldDiff(p)

    Su = p_v .* y_u - p_u .* y_v
    Sv = p_u .* x_v - p_v .* x_u
    # 生成系数
    renew_coff_field!(coodinate, convectionSecheme, n, m, mu, rho, valu, bu, alpha, beta, gamma, Ja, x_u, x_v, y_u, y_v, U, V, u, v, x_uv, y_uv, u, Su, Apu, :u, bounds)

    renew_coff_field!(coodinate, convectionSecheme, n, m, mu, rho, valv, bv, alpha, beta, gamma, Ja, x_u, x_v, y_u, y_v, U, V, u, v, x_uv, y_uv, v, Sv, Apv, :v, bounds)

    lambda=1
    Au = sparse(row, col, valu)
    u = lambda*reshape(Au \ bu, n, m)+(1-lambda)*u
    Av = sparse(row, col, valv)
    v = lambda*reshape(Av \ bv, n, m)+(1-lambda)*v


    valp = Vector{Float64}(undef, 0)
    bp = Vector{Float64}(undef, 0)

    SIMPLE(coodinate, n, m, valp, bp, rho, x_u, x_v, y_u, y_v, Apu, Apv, p, U, V, u, v, x_uv, y_uv, Ja, alpha, gamma, :p, bounds)

    Ap = sparse(rowp, colp, valp)

    dp = reshape(Ap \ bp, n, m)*1.0

    
    dpu, dpv = fieldDiff(dp)
    #u += (-y_v .* dpu + y_u .* dpv) ./ Apu
    #v += (x_v .* dpu - x_u .* dpv) ./ Apv
    
    p += dp
    #=
    ux,uy=fieldDiff(u)
    vx,vy=fieldDiff(v)
    ux+vy
    Ux,Uy=fieldDiff(U)
    Vx,Vy=fieldDiff(V)
    Ux+Vy
    =#
    
    dU=-(y_v.^2 ./ Apu + x_v.^2 ./ Apv).*dpu
    dV=-(y_u.^2 ./ Apu + x_u.^2 ./ Apv).*dpv
    U+=dU
    V+=dV
    Ux,Uy=fieldDiff(U)
    Vx,Vy=fieldDiff(V)
    Ux+Vy
    showFlow(x_uv,y_uv,u,v,0.5)
end

