using EasyCFD2D

l_d=2.0
D_max=15
L_cyc=D_max*l_d
theta=15/180*pi
k_bound=tan(theta)
D_min=D_max-L_cyc*k_bound
bound_control=0.8

function bd1(t)
    [t,0]
end
t1=[0,L_cyc]

#bd2 t in [-1,1]
function bd2(t)
    [L_cyc,D_min/2*(t/D_min*2)^bound_control]
end
t2=[0,D_min/2]

#bd3 t in [0,pi/2]
function bd3(t)
    x=L_cyc-t
    if t<=L_cyc/2
        y=D_min/2+k_bound*t
    else
        y=D_min/2+k_bound*(L_cyc-t)
    end
    return [x,y]
end
t3=[0,L_cyc]

#bd4 t in [2,-2]
function bd4(t)
    [0,D_min/2*(t/D_min*2)^bound_control]
end
t4=[D_min/2,0]

bounds=Vector{bound}(undef,4)
bounds[1]=bound(bd1,t1)
bounds[2]=bound(bd2,t2)
bounds[3]=bound(bd3,t3)
bounds[4]=bound(bd4,t4)

include("PostProcess.jl")

m=50
n=40
@time x_uv,y_uv=GSGrider(m,n,bounds)
gridPlot(x_uv,y_uv)

@time x_uv,y_uv=EasyCFD2D.JacobianGrider(m,n,bounds,maxep=1e-3,relax=0.2,displayStep=1)
gridPlot(x_uv,y_uv)
