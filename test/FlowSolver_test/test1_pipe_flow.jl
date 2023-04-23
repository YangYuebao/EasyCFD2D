using EasyCFD2D
using SparseArrays

function bd1(t)
    [0,1-t^1.5]
end
t1=[0,1]

#bd2 t in [-1,1]
function bd2(t)
    [t,0]
end
t2=[0,5]

#bd3 t in [0,pi/2]
function bd3(t)
    [5,1-(1-t)^1.5]
end
t3=[0,1]

#bd4 t in [2,-2]
function bd4(t)
    [5-t,1]
end
t4=[0,5]

bounds=Vector{bound}(undef,4)
bounds[1]=bound(bd1,t1,pressureInlet(),10)
bounds[2]=bound(bd2,t2,symetryAxis())
bounds[3]=bound(bd3,t3,FDOutlet(),0)
bounds[4]=bound(bd4,t4,stillWall())

include("../Grider_test/PostProcess.jl")

m=5
n=5

x_uv,y_uv=EasyCFD2D.JacobianGrider(m,n,bounds,maxep=1e-3,relax=0.2,displayStep=10)
#gridPlot(x_uv,y_uv)

# 初始化
begin
    u,v,p,S,x_u,y_u,x_v,y_v,alpha,beta,gamma,Ja=EasyCFD2D.fieldA(x_uv,y_uv)
    row,col=EasyCFD2D.index_generation(n,m,EasyCFD2D.getPlace(EasyCFD2D.SecondOrderUpwind()))
    rowp,colp=EasyCFD2D.index_generation(n,m,EasyCFD2D.getPlace(EasyCFD2D.Point5()))

    mu=10.0
    rho=1.1
end

begin    
    val=Vector{Float64}(undef,0)
    b=Vector{Float64}(undef,0)
    Apu=zeros(n,m)
    Apv=zeros(n,m)
end

U,V=EasyCFD2D.getUV(u,v,x_u,x_v,y_u,y_v)
p_u,p_v=EasyCFD2D.fieldDiff(p)

S=p_u.*y_v-p_v.*y_u
#Sv=p_v*x_u-p_u*x_v

# 生成系数
EasyCFD2D.renew_coff_field!(EasyCFD2D.SecondOrderUpwind(), n, m, mu, rho, val, b,alpha,beta,gamma,Ja, U, V, u, S,Apu, :u, bounds)

Apv=Apu[:,:]

A=sparse(row,col,val)
u+=reshape(A\b,n,m)

begin
    valp=Vector{Float64}(undef,0)
    bp=Vector{Float64}(undef,0)
end

EasyCFD2D.SIMPLE(n,m,valp,bp,rho,x_v,y_u,y_v,x_u,Apu,Apv,p,U,V,alpha,gamma,:p,bounds)

B=sparse(rowp,colp,valp)

p+=reshape(B\bp,n,m)

