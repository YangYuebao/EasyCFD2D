abstract type ConvectionScheme <: AbstractScheme end

struct SecondOrderUpwind <:ConvectionScheme end

function getPlace(::SecondOrderUpwind)
    place=Vector{Vector{Int64}}(undef,0)
    for j=-2:2
        for i=-(2-abs(j)):2-abs(j)
            push!(place,[i,j])
        end
    end
    return place
end

