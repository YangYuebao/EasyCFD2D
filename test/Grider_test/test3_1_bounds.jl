using EasyCFD2D
begin
    function bd1(t)
        [3*(t/3)^0.8,0]
    end
    t1=[0,3]

    #bd2 t in [-1,1]
    function bd2(t)
        [3,2*(t*t*(3-2*t))^0.95]
    end
    t2=[0,1]

    #bd3 t in [0,pi/2]
    function bd3(t)
        if t>=0 && t<1.5
            return [3-1.5*(t/1.5)^0.8,2]
        elseif t>=1.5 && t<=2.5
            k=1.1
            x0=0.5
            x1=2
            b=1
            if t<=2.0
                return [1.5,2-(1+1/k)*b/2/x0*(t-x1+k/(k+1)*x0+x0/(k+1)*((x1-t)/x0)^(k+1))]
            else
                return [1.5,1+(1+1/k)*b/2/x0*((4-t)-x1+k/(k+1)*x0+x0/(k+1)*((x1-(4-t))/x0)^(k+1))]
            end
            return [1.5,2-(t-2)]
        elseif t>2.5 && t<=4
            return [1.5-((t-2.5)/1.5)^1*1.5,1]
        end
    end
    t3=[0,4]

    #bd4 t in [2,-2]
    function bd4(t)
        [0,1-t^1.5]
    end
    t4=[0,1]

    bounds=Vector{bound}(undef,4)
    bounds[1]=bound(bd1,t1)
    bounds[2]=bound(bd2,t2)
    bounds[3]=bound(bd3,t3)
    bounds[4]=bound(bd4,t4)
end
include("PostProcess.jl")

m=50
n=40
#@time x_uv,y_uv=GSGrider(m,n,bounds)
#gridPlot(x_uv,y_uv)

@time x_uv,y_uv=EasyCFD2D.JacobianGrider(m,n,bounds,maxep=1e-3,relax=0.2,displayStep=1)
gridPlot(x_uv,y_uv)
