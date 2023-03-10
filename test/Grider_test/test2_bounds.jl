using EasyCFD2D

#bd1 t in [pi/2,0]
function bd1(t)
    [t^3,t]
end
t1=[0,1]

#bd2 t in [-1,1]
function bd2(t)
    [t,2-t]
end
t2=[1,0]

#bd3 t in [0,pi/2]
function bd3(t)
    [-sin(t),1+cos(t)]
end
t3=[0,pi/2]

#bd4 t in [2,-2]
function bd4(t)
    [-sin(t),1+cos(t)]
end
t4=[pi/2,pi]

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

