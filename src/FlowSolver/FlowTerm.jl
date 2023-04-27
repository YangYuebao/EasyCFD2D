# 对流格式
abstract type ConvectionScheme <: AbstractScheme end

struct SecondOrderUpwind <:ConvectionScheme end #二阶迎风

function getPlace(::SecondOrderUpwind)
    place=Vector{Vector{Int64}}(undef,0)
    for j=-2:2
        for i=-(2-abs(j)):2-abs(j)
            push!(place,[i,j])
        end
    end
    return place
end

function index_generation(n::Int64, m::Int64, ::SecondOrderUpwind)
    place = Vector{Vector{Int64}}(undef, 0)
    for j = -2:2
        for i = -(2 - abs(j)):2-abs(j)
            push!(place, [i, j])
        end
    end
    index_generation(n, m, place)
end
