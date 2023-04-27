using EasyCFD2D

begin
    l_d = 2.0
    D_max = 15
    L_cyc = D_max * l_d
    theta = 20 / 180 * pi
    k_bound = tan(theta)
    D_min = D_max - L_cyc * k_bound
    bound_control = 0.8
end
begin
    function bd1(t)
        [0, D_min / 2 * (t / D_min * 2)^bound_control]
    end
    t1 = [D_min / 2, 0]

    function bd2(t)
        [t, 0]
    end
    t2 = [0, L_cyc]

    #bd2 t in [-1,1]
    function bd3(t)
        [L_cyc, D_min / 2 * (t / D_min * 2)^bound_control]
    end
    t3 = [0, D_min / 2]

    #bd3 t in [0,pi/2]
    function bd4(t)
        x = L_cyc - t
        if t <= L_cyc / 2
            y = D_min / 2 + k_bound * t
        else
            y = D_min / 2 + k_bound * x
        end
        return [x, y]
    end
    t4 = [0, L_cyc]

    #bd4 t in [2,-2]

    bounds[1] = bound(bd1, t1, pressureInlet(), 10)
    bounds[2] = bound(bd2, t2, symetryAxis())
    bounds[3] = bound(bd3, t3, FDOutlet(), 0)
    bounds[4] = bound(bd4, t4, stillWall())
end
include("../Grider_test/PostProcess.jl")

m = 41
n = 41
#x_uv,y_uv=GSGrider(m,n,bounds)
#x_uv, y_uv = EasyCFD2D.JacobianGrider(m, n, bounds, maxep=1e-5, relax=0.2, displayStep=1,maxcount=100)

x_uv=zeros(n,m)
y_uv=zeros(n,m)
dx=L_cyc/(m-1)

for j=1:m
    x=dx*(j-1)
    l=bd4(L_cyc-x)[2]
    for i=1:n
        x_uv[i,j]=x
        y_uv[i,j]=l * ( 1 - (i-1) / (n-1) )^bound_control
    end
end

display(gridPlot(x_uv, y_uv))

mu = 100.0
rho = 1.0
uc, vc, pc = solvefield(Rectangular(), SecondOrderUpwind(), x_uv, y_uv, mu, rho, bounds; abstol=1e-5, maxiter=1000)
lambda=1.0
showFlow(x_uv,y_uv,lambda*uc,lambda*vc)