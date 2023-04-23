abstract type AbstractBound end
abstract type AbstractBoundType end

struct stillWall <: AbstractBoundType end
struct movingWall <: AbstractBoundType end
struct FDOutlet <: AbstractBoundType end #fully developed outlet
struct velocityInlet <: AbstractBoundType end
struct pressureInlet <: AbstractBoundType end
struct symetryAxis <: AbstractBoundType end
struct generalBound <: AbstractBoundType end

struct bound <: AbstractBound
    fun::Function
    span::Vector
    u::Function
    v::Function
    p::Function
    bound(fun::Function,span::Vector,u::Function,v::Function,p::Function)=new(fun,span,u,v,p)
end

# 网格划分用例
function bound(fun::Function,span::Vector)
    u(t)=[1,0,0]
    v(t)=[1,0,0]
    p(t)=[1,0,0]
    bound(fun,span,u,v,p)
end

# 静止墙壁边界
function bound(fun::Function,span::Vector,::stillWall)
    u(t)=[1,0,0]
    v(t)=[1,0,0]
    p(t)=[0,1,0]
    bound(fun,span,u,v,p)
end

# 移动墙壁边界
function bound(fun::Function,span::Vector,::movingWall,u0::Number,v0::Number)
    u(t)=[1,0,u0]
    v(t)=[1,0,v0]
    p(t)=[0,1,0]
    bound(fun,span,u,v,p)
end

# 充分发展出口，压力已知
function bound(fun::Function,span::Vector,::FDOutlet,p0::Number)
    u(t)=[0,1,0]
    v(t)=[0,1,0]
    p(t)=[1,0,p0]
    bound(fun,span,u,v,p)
end

# 速度入口，形式上和移动墙壁相同
@deprecate bound(fun::Function,span::Vector,::velocityInlet,u0::Number,v0::Number) bound(fun,span,movingWall(),u0,v0)

# 压力入口，形式上和充分发展出口相同
@deprecate bound(fun::Function,span::Vector,::pressureInlet,p0::Number) bound(fun,span,FDOutlet(),p0)

# 对称轴
function bound(fun::Function,span::Vector,::symetryAxis)
    u(t)=[0,1,0]
    v(t)=[1,0,0]
    p(t)=[0,1,0]
    bound(fun,span,u,v,p)
end

@deprecate bound(fun::Function,span::Vector,::generalBound,u::Function,v::Function,p::Function) bound(fun::Function,span::Vector,u::Function,v::Function,p::Function)

#=
struct extraBounds <: AbstractBound
    name::Symbol
end
=#
