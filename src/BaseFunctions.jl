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

function bound(fun::Function,span::Vector)
    u(t)=[1,0,0]
    v(t)=[1,0,0]
    p(t)=[1,0,0]
    bound(fun,span,u,v,p)
end

struct extraBounds <: AbstractBound
    name::Symbol
end