# 压力速度耦合算法

# SIMPLE算法
function SIMPLE(::Rectangular,n::Int64,m::Int64,val::Vector{Float64},b::Vector{Float64},rho::Float64,x_u::Matrix{Float64},x_v::Matrix{Float64},y_u::Matrix{Float64},y_v::Matrix{Float64},Apu::Matrix{Float64},Apv::Matrix{Float64},p::Matrix{Float64},U::Matrix{Float64},V::Matrix{Float64},u::Matrix{Float64},v::Matrix{Float64},x_uv::Matrix{Float64},y_uv::Matrix{Float64},Ja::Matrix{Float64},alpha::Matrix{Float64},gamma::Matrix{Float64},phi_symbol::Symbol,bounds::Vector{bound})
    for j=1:m
        for i=1:n
            #内点
            if i >= 2 && i <= n - 1 && j >= 2 && j <= m - 1
                val_push = [
                    rho*0.5*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j-1]^2/Apu[i,j-1]+x_v[i,j-1]^2/Apv[i,j-1]),
                    rho*0.5*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i-1,j]^2/Apu[i-1,j]+x_u[i-1,j]^2/Apv[i-1,j]),
                    0,
                    rho*0.5*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i+1,j]^2/Apu[i+1,j]+x_u[i+1,j]^2/Apv[i+1,j]),
                    rho*0.5*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j+1]^2/Apu[i,j+1]+x_v[i,j+1]^2/Apv[i,j+1])
                ]
                val_push[3]=-sum(val_push)
                b_push=rho*0.5*(U[i,j+1]-U[i,j-1]+V[i-1,j]-V[i+1,j])
                push!(val,val_push...)
                push!(b,b_push)
                continue
            end

            # 边界点
            to_val = to_val_index(n,m,i, j,Point5())
            val_push = zeros(5)
            bd_coff = toABC(n,m,i, j,alpha,gamma,phi_symbol,bounds)
            if bd_coff[2] == 0
                val_push[3] = 1
                val_push = val_push[findall(x -> x == true, to_val)]
                push!(val, val_push...)
                push!(b, bd_coff[3]/bd_coff[1]-p[i,j])
                continue
            end
            if i == 1 && j == 1
                val_push = [
                    0,
                    rho*0.25*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i+1,j]^2/Apu[i+1,j]+x_u[i+1,j]^2/Apv[i+1,j]),
                    rho*0.25*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j+1]^2/Apu[i,j+1]+x_v[i,j+1]^2/Apv[i,j+1])
                ]
                val_push[1] = -sum(val_push)
                b_push=0.25*rho*(U[i,j+1]-U[i,j]+V[i,j]-V[i+1,j])
                push!(val, val_push...)
                push!(b, b_push)
                continue
            end
            if i == n && j == 1
                val_push = [
                    rho*0.25*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i-1,j]^2/Apu[i-1,j]+x_u[i-1,j]^2/Apv[i-1,j]),
                    0,
                    rho*0.25*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j+1]^2/Apu[i,j+1]+x_v[i,j+1]^2/Apv[i,j+1])
                ]
                val_push[2] = -sum(val_push)
                b_push=0.25*rho*(U[i,j+1]-U[i,j]+V[i-1,j]-V[i,j])
                push!(val, val_push...)
                push!(b, b_push)
                continue
            end
            if i == n && j == m
                val_push = [
                    rho*0.25*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j-1]^2/Apu[i,j-1]+x_v[i,j-1]^2/Apv[i,j-1]),
                    rho*0.25*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i-1,j]^2/Apu[i-1,j]+x_u[i-1,j]^2/Apv[i-1,j]),
                    0
                ]
                val_push[3] = -sum(val_push)
                b_push=0.25*rho*(U[i,j]-U[i,j-1]+V[i-1,j]-V[i,j])
                push!(val, val_push...)
                push!(b, b_push)
                continue                
            end
            if i == 1 && j == m
                val_push = [
                    rho*0.5*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j-1]^2/Apu[i,j-1]+x_v[i,j-1]^2/Apv[i,j-1]),
                    0,
                    rho*0.5*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i+1,j]^2/Apu[i+1,j]+x_u[i+1,j]^2/Apv[i+1,j]),
                ]
                val_push[2] = -sum(val_push)
                b_push=0.25*rho*(U[i,j]-U[i,j-1]+V[i,j]-V[i+1,j])
                push!(val, val_push...)
                push!(b, b_push)
                continue
            end
            if i == 1
                val_push = [
                    0,
                    1,
                    -1,
                    0
                ]
                push!(val, val_push...)
                push!(b, 0)
                continue
            end
            if i == n
                val_push = [
                    0,
                    -1,
                    1,
                    0
                ]
                push!(val, val_push...)
                push!(b, 0)
                continue
            end
            if j == 1
                val_push = [
                    0,
                    1,
                    0,
                    -1
                ]
                push!(val, val_push...)
                push!(b, 0)
                continue
            end
            if j == m
                val_push = [
                    -1,
                    0,
                    1,
                    0
                ]
                push!(val, val_push...)
                push!(b, 0)
                continue
            end
        end
    end
end

function SIMPLE(::Cylindrical,n::Int64,m::Int64,val::Vector{Float64},b::Vector{Float64},rho::Float64,x_u::Matrix{Float64},x_v::Matrix{Float64},y_u::Matrix{Float64},y_v::Matrix{Float64},Apu::Matrix{Float64},Apv::Matrix{Float64},p::Matrix{Float64},U::Matrix{Float64},V::Matrix{Float64},u::Matrix{Float64},v::Matrix{Float64},x_uv::Matrix{Float64},y_uv::Matrix{Float64},Ja::Matrix{Float64},alpha::Matrix{Float64},gamma::Matrix{Float64},phi_symbol::Symbol,bounds::Vector{bound})
    for j=1:m
        for i=1:n
            S = i==n ? 0 : rho*v[i,j]*Ja[i,j]/y_uv[i,j] # 由坐标线弯曲引起的质量源项
            #内点
            if i >= 2 && i <= n - 1 && j >= 2 && j <= m - 1
                val_push = [
                    rho*0.5*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j-1]^2/Apu[i,j-1]+x_v[i,j-1]^2/Apv[i,j-1]),
                    rho*0.5*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i-1,j]^2/Apu[i-1,j]+x_u[i-1,j]^2/Apv[i-1,j]),
                    0,
                    rho*0.5*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i+1,j]^2/Apu[i+1,j]+x_u[i+1,j]^2/Apv[i+1,j]),
                    rho*0.5*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j+1]^2/Apu[i,j+1]+x_v[i,j+1]^2/Apv[i,j+1])
                ]
                val_push[3]=-sum(val_push)
                b_push=rho*0.5*(U[i,j+1]-U[i,j-1]+V[i-1,j]-V[i+1,j])+S
                push!(val,val_push...)
                push!(b,b_push)
                continue
            end

            # 边界点
            to_val = to_val_index(n,m,i, j,Point5())
            val_push = zeros(5)
            bd_coff = toABC(n,m,i, j,alpha,gamma,phi_symbol,bounds)
            if bd_coff[2] == 0
                val_push[3] = 1
                val_push = val_push[findall(x -> x == true, to_val)]
                push!(val, val_push...)
                push!(b, bd_coff[3]/bd_coff[1]-p[i,j])
                continue
            end
            if i == 1 && j == 1
                val_push = [
                    0,
                    rho*0.25*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i+1,j]^2/Apu[i+1,j]+x_u[i+1,j]^2/Apv[i+1,j]),
                    rho*0.25*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j+1]^2/Apu[i,j+1]+x_v[i,j+1]^2/Apv[i,j+1])
                ]
                val_push[1] = -sum(val_push)
                b_push=0.25*rho*(U[i,j+1]-U[i,j]+V[i,j]-V[i+1,j])+0.25*S
                push!(val, val_push...)
                push!(b, b_push)
                continue
            end
            if i == n && j == 1
                val_push = [
                    rho*0.25*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i-1,j]^2/Apu[i-1,j]+x_u[i-1,j]^2/Apv[i-1,j]),
                    0,
                    rho*0.25*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j+1]^2/Apu[i,j+1]+x_v[i,j+1]^2/Apv[i,j+1])
                ]
                val_push[2] = -sum(val_push)
                b_push=0.25*rho*(U[i,j+1]-U[i,j]+V[i-1,j]-V[i,j])+0.25*S
                push!(val, val_push...)
                push!(b, b_push)
                continue
            end
            if i == n && j == m
                val_push = [
                    rho*0.25*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j-1]^2/Apu[i,j-1]+x_v[i,j-1]^2/Apv[i,j-1]),
                    rho*0.25*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i-1,j]^2/Apu[i-1,j]+x_u[i-1,j]^2/Apv[i-1,j]),
                    0
                ]
                val_push[3] = -sum(val_push)
                b_push=0.25*rho*(U[i,j]-U[i,j-1]+V[i-1,j]-V[i,j])+0.25*S
                push!(val, val_push...)
                push!(b, b_push)
                continue                
            end
            if i == 1 && j == m
                val_push = [
                    rho*0.5*(y_v[i,j]^2/Apu[i,j]+x_v[i,j]^2/Apv[i,j]+y_v[i,j-1]^2/Apu[i,j-1]+x_v[i,j-1]^2/Apv[i,j-1]),
                    0,
                    rho*0.5*(y_u[i,j]^2/Apu[i,j]+x_u[i,j]^2/Apu[i,j]+y_u[i+1,j]^2/Apu[i+1,j]+x_u[i+1,j]^2/Apv[i+1,j]),
                ]
                val_push[2] = -sum(val_push)
                b_push=0.25*rho*(U[i,j]-U[i,j-1]+V[i,j]-V[i+1,j])+0.25*S
                push!(val, val_push...)
                push!(b, b_push)
                continue
            end
            if i == 1
                val_push = [
                    0,
                    1,
                    -1,
                    0
                ]
                push!(val, val_push...)
                push!(b, 0.5*S)
                continue
            end
            if i == n
                val_push = [
                    0,
                    -1,
                    1,
                    0
                ]
                push!(val, val_push...)
                push!(b, 0.5*S)
                continue
            end
            if j == 1
                val_push = [
                    0,
                    1,
                    0,
                    -1
                ]
                push!(val, val_push...)
                push!(b, 0.5*S)
                continue
            end
            if j == m
                val_push = [
                    -1,
                    0,
                    1,
                    0
                ]
                push!(val, val_push...)
                push!(b, 0.5*S)
                continue
            end
        end
    end
end