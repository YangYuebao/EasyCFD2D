function index_generation(n::Int64, m::Int64, place::Vector{Vector{Int64}})
    row = Vector{Int64}(undef, 0) # 不变
    col = Vector{Int64}(undef, 0) # 不变
    for j = 1:m
        for i = 1:n
            for move in place
                temp = [i, j] + move
                if temp[1] >= 1 && temp[1] <= n && temp[2] >= 1 && temp[2] <= m
                    push!(row, ij2k(n, m, i, j))
                    push!(col, ij2k(n, m, temp[1], temp[2]))
                end
            end
        end
    end

    return row, col
end

function k2ij(n::Int64, m::Int64, k::Int64)
    j::Int64 = 1 + floor(Int64, k / n)
    i::Int64 = k % n
    return [i, j]
end

function ij2k(n::Int64, m::Int64, i::Int64, j::Int64)
    i + (j - 1) * n
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

function renew_coff_field!(::SecondOrderUpwind, n::Int64, m::Int64,mu::Float64,rho::Float64, val::Vector{Float64}, b::Vector{Float64}, U::Vector{Float64}, V::Vector{Float64}, phi::Vector{Float64},phi_symbol::Symbol,bounds::Vector{bound})
    for j = 1:m
        for i = 1:n
            if i >= 3 && i <= n - 2 && j >= 3 && j <= m - 2 #内点
                val_push = -[
                    0,
                    0.125 * mu * (2 * beta[i, j] / Ja[i, j] + beta[i, j-1] / Ja[i, j-1] + beta[i-1, j] / Ja[i-1, j]),
                    0.125 * mu * (beta[i-1, j] / Ja[i-1, j] - beta[i+1, j] / Ja[i+1, j]),
                    -0.125 * mu * (2 * beta[i, j] / Ja[i, j] + beta[i+1, j] / Ja[i+1, j] + beta[i, j-1] / Ja[i, j-1]),
                    0,
                    -0.125 * mu * (beta[i, j+1] / Ja[i, j+1] - beta[i, j-1] / Ja[i, j-1]),
                    0,
                    0.125 * mu * (beta[i, j+1] / Ja[i, j+1] - beta[i, j-1] / Ja[i, j-1]),
                    0,
                    -0.125 * mu * (2 * beta[i, j] / Ja[i, j] + beta[i, j+1] / Ja[i, j+1] + beta[i-1, j] / Ja[i-1, j]),
                    -0.125 * mu * (beta[i-1, j] / Ja[i-1, j] - beta[i+1, j] / Ja[i+1, j]),
                    0.125 * mu * (2 * beta[i, j] / Ja[i, j] + beta[i, j+1] / Ja[i, j+1] + beta[i+1, j] / Ja[i+1, j]),
                    0
                ]
                val_push[3] -= 0.5 * mu * (alpha[i, j] / Ja[i, j] + alpha[i, j-1] / Ja[i, j-1])
                val_push[6] -= 0.5 * mu * (alpha[i, j] / Ja[i, j] + alpha[i-1, j] / Ja[i-1, j])
                val_push[8] -= 0.5 * mu * (alpha[i, j] / Ja[i, j] + alpha[i+1, j] / Ja[i+1, j])
                val_push[11] -= 0.5 * mu * (alpha[i, j] / Ja[i, j] + alpha[i, j+1] / Ja[i, j+1])
                val_push[7] = val_push[3] + val_push[6] + val_push[8] + val_push[11]

                b_push = S[i, j]
                if U[i, j] + U[i, j+1] > 0 && U[i, j-1] + U[i, j] > 0
                    val_push[1] += rho * U[i, j] * (0.5)
                    val_push[3] += rho * U[i, j] * (-2)
                    val_push[7] += rho * U[i, j] * (1.5)
                elseif U[i, j] + U[i, j+1] < 0 && U[i, j-1] + U[i, j] < 0
                    val_push[7] += rho * U[i, j] * (-1.5)
                    val_push[11] += rho * U[i, j] * (2)
                    val_push[13] += rho * U[i, j] * (-0.5)
                elseif U[i, j] + U[i, j+1] >= 0 && U[i, j-1] + U[i, j] <= 0
                    val_push[3] += rho * U[i, j] * (-0.5)
                    val_push[11] += rho * U[i, j] * (0.5)
                else
                    val_push[1] += rho * U[i, j] * (0.5)
                    val_push[3] += rho * U[i, j] * (-1.5)
                    val_push[11] += rho * U[i, j] * (1.5)
                    val_push[13] += rho * U[i, j] * (-0.5)
                end
                if V[i, j] + V[i-1, j] > 0 && V[i, j] + V[i+1, j] > 0
                    val_push[7] += rho * V[i, j] * (1.5)
                    val_push[8] += rho * V[i, j] * (-2)
                    val_push[9] += rho * V[i, j] * (0.5)
                elseif V[i, j] + V[i-1, j] < 0 && V[i, j] + V[i-1, j] < 0
                    val_push[5] += rho * V[i, j] * (-0.5)
                    val_push[6] += rho * V[i, j] * (2)
                    val_push[7] += rho * V[i, j] * (-1.5)
                elseif V[i, j] + V[i, j+1] <= 0 && V[i, j-1] + V[i, j] >= 0
                    val_push[6] += rho * V[i, j] * (+0.5)
                    val_push[8] += rho * V[i, j] * (-0.5)
                else
                    val_push[5] += rho * V[i, j] * (-0.5)
                    val_push[6] += rho * V[i, j] * (+1.5)
                    val_push[8] += rho * V[i, j] * (-1.5)
                    val_push[9] += rho * V[i, j] * (+0.5)
                end
                push!(val, val_push)
                push!(b, b_push)
                continue
            end

            ## 处理边界节点
            to_val = to_val_index(i, j)
            val_push = zeros(13)
            if (i == 2 || j == 2 || i == n - 1 || j == m - 1)&&(i!=1)&&(i!=n)&&(j!=1)&&(j!=m) #二层外点
                val_push = -[
                    0,
                    0.125 * mu * (2 * beta[i, j] / Ja[i, j] + beta[i, j-1] / Ja[i, j-1] + beta[i-1, j] / Ja[i-1, j]),
                    0.125 * mu * (beta[i-1, j] / Ja[i-1, j] - beta[i+1, j] / Ja[i+1, j]),
                    -0.125 * mu * (2 * beta[i, j] / Ja[i, j] + beta[i+1, j] / Ja[i+1, j] + beta[i, j-1] / Ja[i, j-1]),
                    0,
                    -0.125 * mu * (beta[i, j+1] / Ja[i, j+1] - beta[i, j-1] / Ja[i, j-1]),
                    0,
                    0.125 * mu * (beta[i, j+1] / Ja[i, j+1] - beta[i, j-1] / Ja[i, j-1]),
                    0,
                    -0.125 * mu * (2 * beta[i, j] / Ja[i, j] + beta[i, j+1] / Ja[i, j+1] + beta[i-1, j] / Ja[i-1, j]),
                    -0.125 * mu * (beta[i-1, j] / Ja[i-1, j] - beta[i+1, j] / Ja[i+1, j]),
                    0.125 * mu * (2 * beta[i, j] / Ja[i, j] + beta[i, j+1] / Ja[i, j+1] + beta[i+1, j] / Ja[i+1, j]),
                    0
                ]
                val_push[3] -= 0.5 * mu * (alpha[i, j] / Ja[i, j] + alpha[i, j-1] / Ja[i, j-1])
                val_push[6] -= 0.5 * mu * (alpha[i, j] / Ja[i, j] + alpha[i-1, j] / Ja[i-1, j])
                val_push[8] -= 0.5 * mu * (alpha[i, j] / Ja[i, j] + alpha[i+1, j] / Ja[i+1, j])
                val_push[11] -= 0.5 * mu * (alpha[i, j] / Ja[i, j] + alpha[i, j+1] / Ja[i, j+1])
                val_push[7] = val_push[3] + val_push[6] + val_push[8] + val_push[11]

                b_push = S[i, j]

                if U[i, j] + U[i, j+1] > 0 && U[i, j-1] + U[i, j] > 0
                    if to_val[1]==true
                        val_push[1] += rho * U[i, j] * (0.5)
                    else
                        b_push -= rho * U[i, j] * (0.5) * (2*phi[i,1]-phi[i,2])
                    end
                    val_push[3] += rho * U[i, j] * (-2)
                    val_push[7] += rho * U[i, j] * (1.5)
                elseif U[i, j] + U[i, j+1] < 0 && U[i, j-1] + U[i, j] < 0
                    val_push[7] += rho * U[i, j] * (-1.5)
                    val_push[11] += rho * U[i, j] * (2)
                    if to_val[13]==true
                        val_push[13] += rho * U[i, j] * (-0.5)
                    else
                        b_push -= rho * U[i, j] * (-0.5) * (2*phi[i,m]-phi[i,m-1])
                    end
                elseif U[i, j] + U[i, j+1] >= 0 && U[i, j-1] + U[i, j] <= 0
                    val_push[3] += rho * U[i, j] * (-0.5)
                    val_push[11] += rho * U[i, j] * (0.5)
                else
                    if to_val[1]==true
                        val_push[1] += rho * U[i, j] * (0.5)
                    else
                        b_push -= rho * U[i, j] * (0.5) * (2*phi[i,1]-phi[i,2])
                    end 
                    val_push[3] += rho * U[i, j] * (-1.5)
                    val_push[11] += rho * U[i, j] * (1.5)
                    if to_val[13]==true
                        val_push[13] += rho * U[i, j] * (-0.5)
                    else
                        b_push -= rho * U[i, j] * (-0.5) * (2*phi[i,m]-phi[i,m-1])
                    end
                end
                if V[i, j] + V[i-1, j] > 0 && V[i, j] + V[i+1, j] > 0
                    val_push[7] += rho * V[i, j] * (1.5)
                    val_push[8] += rho * V[i, j] * (-2)
                    if to_val[9]==true
                        val_push[9] += rho * V[i, j] * (0.5)
                    else
                        b_push -= rho * V[i, j] * (0.5)*(2*phi[n,j]-phi[n-1,j])
                    end
                elseif V[i, j] + V[i-1, j] < 0 && V[i, j] + V[i-1, j] < 0
                    if to_val[5]==true
                        val_push[5] += rho * V[i, j] * (-0.5)
                    else
                        b_push -= rho * V[i, j] * (-0.5)*(2*phi[n,j]-phi[n-1,j])
                    end
                    val_push[6] += rho * V[i, j] * (2)
                    val_push[7] += rho * V[i, j] * (-1.5)
                elseif V[i, j] + V[i, j+1] <= 0 && V[i, j-1] + V[i, j] >= 0
                    val_push[6] += rho * V[i, j] * (+0.5)
                    val_push[8] += rho * V[i, j] * (-0.5)
                else
                    if to_val[5]==true
                        val_push[5] += rho * V[i, j] * (-0.5)
                    else
                        b_push -= rho * V[i, j] * (-0.5)*(2*phi[n,j]-phi[n-1,j])
                    end
                    val_push[6] += rho * V[i, j] * (+1.5)
                    val_push[8] += rho * V[i, j] * (-1.5)
                    if to_val[9]==true
                        val_push[9] += rho * V[i, j] * (0.5)
                    else
                        b_push -= rho * V[i, j] * (0.5)*(2*phi[n,j]-phi[n-1,j])
                    end
                end
                val_push = val_push[findall(x->x==true,res)]
                push!(val,val_push)
                push!(b,b_push)
                continue
            end


            ## 处理最外层边界节点
            
        end
    end

    function to_val_index(i::Int64, j::Int64)
        res = ones(Bool, 13)
        place = getPlace(SecondOrderUpwind())
        for k = 1:13
            place[k] += [i, j]
            if place[k][1] < 1 || place[k][1] > n || place[k][2] < 1 || place[k][2] > m
                res[k] = false
            end
        end
        return res
    end

    function toABC(i::Int64,j::Int64)
        if i==1 && j==1
            
    end
end