function index_generation(n::Int64,m::Int64,place::Vector{Vector{Int64}})    
    row = Vector{Int64}(undef, 0) # 不变
    col = Vector{Int64}(undef, 0) # 不变
    for j=1:m
        for i=1:n        
            for move in place
                temp=[i,j]+move
                if temp[1]>=1 && temp[1]<=n && temp[2]>=1 && temp[2]<=m
                    push!(row,ij2k(n,m,i,j))
                    push!(col,ij2k(n,m,temp[1],temp[2]))
                end
            end
        end
    end

    return row,col
end

function k2ij(n::Int64,m::Int64,k::Int64)
    j::Int64=1+floor(Int64,k/n)
    i::Int64=k % n
    return [i,j]
end

function ij2k(n::Int64,m::Int64,i::Int64,j::Int64)
    i+(j-1)*n
end

function index_generation(n::Int64,m::Int64,::SecondOrderUpwind)
    place=Vector{Vector{Int64}}(undef,0)
    for j=-2:2
        for i=-(2-abs(j)):2-abs(j)
            push!(place,[i,j])
        end
    end
    index_generation(n,m,place)
end

function renew_coff_field!()
    
end