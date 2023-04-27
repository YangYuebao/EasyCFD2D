using EasyCFD2D,SparseArrays

begin
    L=1.0
    mu=1/100
    rho=1.0
    u0=0.01 
end

begin
    function myfun(x,x1,y1)
        if x>=0 && x<=x1
            return y1/x1*x
        elseif x>1-x1 && x<=1
            return 1-y1/x1*(1-x)
        else
            return 0.5+(x-0.5)*(1-2*y1)/(1-2*x1)
        end
    end
    x1=0.2
    y1=0.05
    function bd1(t)
        [0, L*(1 - myfun(t/L,x1,y1))]
    end
    t1 = [0, L]

    #bd2 t in [-1,1]
    function bd2(t)
        [L*myfun(t/L,x1,y1), 0]
    end
    t2 = [0, L]

    #bd3 t in [0,pi/2]
    function bd3(t)
        [L, L*myfun(t/L,x1,y1)]
    end
    t3 = [0, L]

    #bd4 t in [2,-2]
    function bd4(t)
        [L*(1 - myfun(t/L,x1,y1)), L]
    end
    t4 = [0, L]

    bounds = Vector{bound}(undef, 4)
    bounds[1] = bound(bd1, t1, stillWall())
    bounds[2] = bound(bd2, t2, stillWall())
    bounds[3] = bound(bd3, t3, stillWall())
    bounds[4] = bound(bd4, t4, movingWall(),u0,0)
end

begin
    m = 35
    n = 35
    #x_uv,y_uv=GSGrider(m,n,bounds)
    x_uv, y_uv = EasyCFD2D.JacobianGrider(n, m, bounds, maxep=1e-3, relax=0.2, displayStep=10)
    gridPlot(x_uv,y_uv)
end
begin
    coodinate=Rectangular()
    convectionSecheme=SecondOrderUpwind()
    u, v, p, x_u, y_u, x_v, y_v, alpha, beta, gamma, Ja = fieldA(x_uv, y_uv)
    row, col = EasyCFD2D.index_generation(n, m, getPlace(SecondOrderUpwind()))
    rowp, colp = EasyCFD2D.index_generation(n, m, getPlace(EasyCFD2D.Point5()))
end

begin
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

    lambda=1.0
    A = sparse(row, col, valu)
    u = lambda*reshape(A \ bu, n, m)+(1-lambda)*u
    A = sparse(row, col, valv)
    v = lambda*reshape(A \ bv, n, m)+(1-lambda)*v


    valp = Vector{Float64}(undef, 0)
    bp = Vector{Float64}(undef, 0)

    SIMPLE(coodinate, n, m, valp, bp, rho, x_u,x_v,y_u,y_v, Apu, Apv, p, U, V,u,v,x_uv,y_uv,Ja, alpha, gamma, :p, bounds)

    Ap = sparse(rowp, colp, valp)

    dp = reshape(Ap \ bp, n, m)
    #show(Matrix(Ap))
    dpu,dpv=fieldDiff(dp)
    U += -(y_v.^2 ./ Apu + x_v.^2 ./ Apv).*dpu
    V += -(y_u.^2 ./ Apu + x_u.^2 ./ Apv).*dpv
    p+=dp
    
    Ux,Uy=fieldDiff(U)
    Vx,Vy=fieldDiff(V)
    Ux+Vy
    #showFlow(x_uv,y_uv,u,v,0.1)
end

output([Matrix(Ap) bp])
showBound(bounds[4],10)