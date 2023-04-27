using EasyCFD2D

begin
    L=1.0
    mu=1.0
    rho=1.0
    u0=1.0 
end

begin
    function bd1(t)
        [0, L - t]
    end
    t1 = [0, L]

    #bd2 t in [-1,1]
    function bd2(t)
        [t, 0]
    end
    t2 = [0, L]

    #bd3 t in [0,pi/2]
    function bd3(t)
        [L, t]
    end
    t3 = [0, L]

    #bd4 t in [2,-2]
    function bd4(t)
        [L - t, L]
    end
    t4 = [0, L]

    bounds = Vector{bound}(undef, 4)
    bounds[1] = bound(bd1, t1, stillWall())
    bounds[2] = bound(bd2, t2, stillWall())
    bounds[3] = bound(bd3, t3, stillWall())
    bounds[4] = bound(bd4, t4, movingWall(),u0,0)
end

begin
    m = 20
    n = 20
    #x_uv,y_uv=GSGrider(m,n,bounds)
    x_uv, y_uv = EasyCFD2D.JacobianGrider(n, m, bounds, maxep=1e-3, relax=0.2, displayStep=10)
    #gridPlot(x_uv,y_uv)
end

ur, vr, pr = solvefield(Rectangular(), SecondOrderUpwind(), x_uv, y_uv, mu, rho, bounds; abstol=1e-5, maxiter=100)
display(quiver(x_uv,y_uv,quiver=(ur,vr)))
