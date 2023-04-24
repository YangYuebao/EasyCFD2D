using EasyCFD2D
using SparseArrays

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
    bounds[1] = bound(bd1, t1, pressureInlet(), 10)
    bounds[2] = bound(bd2, t2, symetryAxis())
    bounds[3] = bound(bd3, t3, FDOutlet(), 0)
    bounds[4] = bound(bd4, t4, stillWall())
end
#include("../Grider_test/PostProcess.jl")

m = 10
n = 10
#x_uv,y_uv=GSGrider(m,n,bounds)
x_uv, y_uv = EasyCFD2D.JacobianGrider(n, m, bounds, maxep=1e-3, relax=0.2, displayStep=10)
#gridPlot(x_uv,y_uv)

# 初始化
begin
    u, v, p, x_u, y_u, x_v, y_v, alpha, beta, gamma, Ja = EasyCFD2D.fieldA(x_uv, y_uv)
    row, col = EasyCFD2D.index_generation(n, m, EasyCFD2D.getPlace(EasyCFD2D.SecondOrderUpwind()))
    rowp, colp = EasyCFD2D.index_generation(n, m, EasyCFD2D.getPlace(EasyCFD2D.Point5()))

    mu = 1.0
    rho = 1.0

    convectionSecheme = SecondOrderUpwind()
    coodinate = Rectangular()#Cylindrical()
    u, v, p, x_u, y_u, x_v, y_v, alpha, beta, gamma, Ja = EasyCFD2D.fieldA(x_uv, y_uv)
    row, col = EasyCFD2D.index_generation(n, m, EasyCFD2D.getPlace(convectionSecheme))
    rowp, colp = EasyCFD2D.index_generation(n, m, EasyCFD2D.getPlace(EasyCFD2D.Point5()))
end


begin
    begin
        valu = Vector{Float64}(undef, 0)
        bu = Vector{Float64}(undef, 0)
        valv = Vector{Float64}(undef, 0)
        bv = Vector{Float64}(undef, 0)
        Apu = zeros(n, m)
        Apv = zeros(n, m)
        U, V = getUV(u, v, x_u, x_v, y_u, y_v)
        p_u, p_v = EasyCFD2D.fieldDiff(p)
    end


    Su = p_v .* y_u - p_u .* y_v
    temp_su = Su[:, :]
    Sv = p_u .* x_v - p_v .* x_u

    #EasyCFD2D.sorceTerm!(coodinate, Su, rho, mu, u, v, y_uv, x_u, x_v, Ja)

    phi_u,phi_v=fieldDiff(u)
    temp = ( - rho * v[1:end-1,:] .* u[1:end-1,:] ./ y_uv[1:end-1,:]) .* Ja[1:end-1,:] + mu ./ y_uv[1:end-1,:] .* (-x_v[1:end-1,:] .* phi_u[1:end-1,:] + x_u[1:end-1,:] .* phi_v[1:end-1,:])
    Su[1:end-1,:] += temp
    Su[end,:] += 2*temp[end,:] - temp[end-1,:]

    EasyCFD2D.sorceTerm!(coodinate, Sv, rho, mu, v, v, y_uv, x_u, x_v, Ja)

    # 生成系数
    EasyCFD2D.renew_coff_field!(convectionSecheme, n, m, mu, rho, valu, bu, alpha, beta, gamma, Ja, U, V, u, Su, Apu, :u, bounds)

    EasyCFD2D.renew_coff_field!(convectionSecheme, n, m, mu, rho, valv, bv, alpha, beta, gamma, Ja, U, V, v, Sv, Apv, :v, bounds)

    A = sparse(row, col, valu)
    u = reshape(A \ bu, n, m)
    A = sparse(row, col, valv)
    v = reshape(A \ bv, n, m)


    valp = Vector{Float64}(undef, 0)
    bp = Vector{Float64}(undef, 0)
    Sp = zeros(n, m)
    EasyCFD2D.sorceTerm!(coodinate, Sp, rho, mu, ones(n, m), v, y_uv, x_u, x_v, Ja)
    SIMPLE(n, m, valp, bp, rho, x_v, y_u, y_v, x_u, Apu, Apv, p, U, V, Sp, alpha, gamma, :p, bounds)

    A = sparse(rowp, colp, valp)

    p += reshape(A \ bp, n, m)
    u
end

#=
phi_field=u[:,:]
phi_u, phi_v = EasyCFD2D.fieldDiff(phi_field)
temp = (-rho * v[1:end-1, :] .* phi_field[1:end-1, :] ./ y_uv[1:end-1, :]) .* Ja[1:end-1, :] + mu ./ y_uv[1:end-1, :] .* (-x_v[1:end-1, :] .* phi_u[1:end-1, :] + x_u[1:end-1, :] .* phi_v[1:end-1, :])
=#
