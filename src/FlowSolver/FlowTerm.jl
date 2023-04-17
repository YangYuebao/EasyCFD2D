abstract type AbstractConvectionTermScheme end

struct SecondOrderUpwind <:AbstractConvectionTermScheme end

function getPlace(::SecondOrderUpwind)
    place=Vector{Vector{Int64}}(undef,0)
    for j=-2:2
        for i=-(2-abs(j)):2-abs(j)
            push!(place,[i,j])
        end
    end
    return place
end

function getCoff(::SecondOrderUpwind)
    coff=zeros(13)
    
end
