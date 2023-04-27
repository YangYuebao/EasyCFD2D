#==============================#
using SparseArrays, EasyCFD2D

begin
    con = 1.0
    function bd1(t)
        [0, 1 - t^con]
    end
    t1 = [0, 1]

    #bd2 t in [-1,1]
    function bd2(t)
        [t, 0]
    end
    t2 = [0, 5]

    #bd3 t in [0,pi/2]
    function bd3(t)
        [5, 1 - (1 - t)^con]
    end
    t3 = [0, 1]

    #bd4 t in [2,-2]
    function bd4(t)
        [5 - t, 1]
    end
    t4 = [0, 5]

    bounds = Vector{bound}(undef, 4)
    bounds[1] = bound(bd1, t1, stillWall())
    bounds[2] = bound(bd2, t2, pressureInlet(), 40)
    bounds[3] = bound(bd3, t3, stillWall())
    bounds[4] = bound(bd4, t4, FDOutlet(), 0)
end

include("../Grider_test/PostProcess.jl")

begin
    m = 21
    n = 21
    #x_uv,y_uv=GSGrider(m,n,bounds)
    x_uv, y_uv = EasyCFD2D.JacobianGrider(n, m, bounds, maxep=1e-3, relax=0.2, displayStep=10)
    #gridPlot(x_uv,y_uv)
    mu = 1.0
    rho = 1.0
    convectionSecheme=SecondOrderUpwind()
    coodinate=Rectangular()
end


#uc, vc, pc = solvefield(Cylindrical(), SecondOrderUpwind(), x_uv, y_uv, mu, rho, bounds; abstol=1e-5, maxiter=100)
#display(quiver(x_uv,y_uv,quiver=(uc,vc)))

begin
    u, v, p, x_u, y_u, x_v, y_v, alpha, beta, gamma, Ja = fieldA(x_uv, y_uv)
    row, col = EasyCFD2D.index_generation(n, m, getPlace(convectionSecheme))
    rowp, colp = EasyCFD2D.index_generation(n, m, getPlace(EasyCFD2D.Point5()))

end
#=
    mu = 1.0
    rho = 1.0
=#
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
    renew_coff_field!(coodinate, convectionSecheme, n, m, mu, rho, valu, bu, alpha, beta, gamma, Ja, x_u, x_v, y_u, y_v, U, V, u, v, x_uv, y_uv, u, Su, Apu, :u, bounds)

    renew_coff_field!(coodinate, convectionSecheme, n, m, mu, rho, valv, bv, alpha, beta, gamma, Ja, x_u, x_v, y_u, y_v, U, V, u, v, x_uv, y_uv, v, Sv, Apv, :v, bounds)

    A = sparse(row, col, valu)
    u = reshape(A \ bu, n, m)
    A = sparse(row, col, valv)
    v = reshape(A \ bv, n, m)


    valp = Vector{Float64}(undef, 0)
    bp = Vector{Float64}(undef, 0)

    SIMPLE(coodinate, n, m, valp, bp, rho, x_u, x_v, y_u, y_v, Apu, Apv, p, U, V, u, v, x_uv, y_uv, Ja, alpha, gamma, :p, bounds)

    A = sparse(rowp, colp, valp)

    dp = reshape(A \ bp, n, m)

        
    dpu,dpv=fieldDiff(dp)
    U += -(y_v.^2 ./ Apu + x_v.^2 ./ Apv).*dpu
    V += -(y_u.^2 ./ Apu + x_u.^2 ./ Apv).*dpv
    p+=dp

end
maximum(v)
#showFlow(x_uv, y_uv, ur, vr, 0.002)
