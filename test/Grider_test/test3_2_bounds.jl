using EasyCFD2D

function bd1(t)
    [t,0]
end
t1=[0,3]

#bd2 t in [-1,1]
function bd2(t)
    [3,t]
end
t2=[0,2]

#bd3 t in [0,pi/2]
function bd3(t)
    if t>=0 && t<1.5
        return [3-t,2]
    elseif t>=1.5 && t<=2.5
        return [1.5,2-t+1.5]
    elseif t>2.5 && t<=4
        return [4-t,1]
    end
end
t3=[0,4]

#bd4 t in [2,-2]
function bd4(t)
    [0,1-t]
end
t4=[0,1]

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
