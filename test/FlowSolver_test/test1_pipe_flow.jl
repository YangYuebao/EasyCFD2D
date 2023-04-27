#==============================#
using SparseArrays,EasyCFD2D

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

begin
    m = 20
    n = 20
    #x_uv,y_uv=GSGrider(m,n,bounds)
    x_uv, y_uv = EasyCFD2D.JacobianGrider(n, m, bounds, maxep=1e-3, relax=0.2, displayStep=10)
    #gridPlot(x_uv,y_uv)
end

mu = 1.0
rho = 1.0
uc, vc, pc = solvefield(Cylindrical(), SecondOrderUpwind(), x_uv, y_uv, mu, rho, bounds; abstol=1e-5, maxiter=100)
display(quiver(x_uv,y_uv,quiver=(uc,vc)))

ur, vr, pr = solvefield(Rectangular(), SecondOrderUpwind(), x_uv, y_uv, mu, rho, bounds; abstol=1e-5, maxiter=100)
display(quiver(x_uv,y_uv,quiver=(ur,vr)))
