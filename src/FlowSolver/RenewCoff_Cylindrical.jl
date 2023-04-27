
function renew_coff_field!(::Cylindrical,::SecondOrderUpwind, n::Int64, m::Int64, mu::Float64, rho::Float64, val::Vector{Float64}, b::Vector{Float64},alpha::Matrix{Float64},beta::Matrix{Float64},gamma::Matrix{Float64},Ja::Matrix{Float64},x_u::Matrix{Float64},x_v::Matrix{Float64},y_u::Matrix{Float64},y_v::Matrix{Float64}, U::Matrix{Float64}, V::Matrix{Float64},u::Matrix{Float64}, v::Matrix{Float64}, x_uv::Matrix{Float64}, y_uv::Matrix{Float64}, phi::Matrix{Float64}, S::Matrix{Float64},Ap::Matrix{Float64}, phi_symbol::Symbol, bounds::Vector{bound})
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
                val_push[6] -= 0.5 * mu * (gamma[i, j] / Ja[i, j] + gamma[i-1, j] / Ja[i-1, j])
                val_push[8] -= 0.5 * mu * (gamma[i, j] / Ja[i, j] + gamma[i+1, j] / Ja[i+1, j])
                val_push[11] -= 0.5 * mu * (alpha[i, j] / Ja[i, j] + alpha[i, j+1] / Ja[i, j+1])
                val_push[7] = -val_push[3] - val_push[6] - val_push[8] - val_push[11]
                
                # 柱坐标系的系数不同处
                val_push[7]+=rho*v[i,j]/y_uv[i,j]*Ja[i,j]
                val_push[3]-=mu*x_v[i,j]/y_uv[i,j]/2
                val_push[6]-=mu*x_u[i,j]/y_uv[i,j]/2
                val_push[8]+=mu*x_u[i,j]/y_uv[i,j]/2
                val_push[11]+=mu*x_v[i,j]/y_uv[i,j]/2

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
                Ap[i,j]=val_push[7]
                push!(val, val_push...)
                push!(b, b_push)
                continue
            end

            ## 处理边界节点
            to_val = to_val_index(n,m,i, j, SecondOrderUpwind())
            val_push = zeros(13)
            #二层外点
            if (i == 2 || j == 2 || i == n - 1 || j == m - 1) && (i != 1) && (i != n) && (j != 1) && (j != m)
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
                val_push[6] -= 0.5 * mu * (gamma[i, j] / Ja[i, j] + gamma[i-1, j] / Ja[i-1, j])
                val_push[8] -= 0.5 * mu * (gamma[i, j] / Ja[i, j] + gamma[i+1, j] / Ja[i+1, j])
                val_push[11] -= 0.5 * mu * (alpha[i, j] / Ja[i, j] + alpha[i, j+1] / Ja[i, j+1])
                val_push[7] = -val_push[3] - val_push[6] - val_push[8] - val_push[11]

                val_push[7]+=rho*v[i,j]/y_uv[i,j]*Ja[i,j]
                val_push[3]-=mu*x_v[i,j]/y_uv[i,j]/2
                val_push[6]-=mu*x_u[i,j]/y_uv[i,j]/2
                val_push[8]+=mu*x_u[i,j]/y_uv[i,j]/2
                val_push[11]+=mu*x_v[i,j]/y_uv[i,j]/2

                b_push = S[i, j]

                if U[i, j] + U[i, j+1] > 0 && U[i, j-1] + U[i, j] > 0
                    if to_val[1] == true
                        val_push[1] += rho * U[i, j] * (0.5)
                    else
                        b_push -= rho * U[i, j] * (0.5) * (2 * phi[i, 1] - phi[i, 2])
                    end
                    val_push[3] += rho * U[i, j] * (-2)
                    val_push[7] += rho * U[i, j] * (1.5)
                elseif U[i, j] + U[i, j+1] < 0 && U[i, j-1] + U[i, j] < 0
                    val_push[7] += rho * U[i, j] * (-1.5)
                    val_push[11] += rho * U[i, j] * (2)
                    if to_val[13] == true
                        val_push[13] += rho * U[i, j] * (-0.5)
                    else
                        b_push -= rho * U[i, j] * (-0.5) * (2 * phi[i, m] - phi[i, m-1])
                    end
                elseif U[i, j] + U[i, j+1] >= 0 && U[i, j-1] + U[i, j] <= 0
                    val_push[3] += rho * U[i, j] * (-0.5)
                    val_push[11] += rho * U[i, j] * (0.5)
                else
                    if to_val[1] == true
                        val_push[1] += rho * U[i, j] * (0.5)
                    else
                        b_push -= rho * U[i, j] * (0.5) * (2 * phi[i, 1] - phi[i, 2])
                    end
                    val_push[3] += rho * U[i, j] * (-1.5)
                    val_push[11] += rho * U[i, j] * (1.5)
                    if to_val[13] == true
                        val_push[13] += rho * U[i, j] * (-0.5)
                    else
                        b_push -= rho * U[i, j] * (-0.5) * (2 * phi[i, m] - phi[i, m-1])
                    end
                end
                if V[i, j] + V[i-1, j] > 0 && V[i, j] + V[i+1, j] > 0
                    val_push[7] += rho * V[i, j] * (1.5)
                    val_push[8] += rho * V[i, j] * (-2)
                    if to_val[9] == true
                        val_push[9] += rho * V[i, j] * (0.5)
                    else
                        b_push -= rho * V[i, j] * (0.5) * (2 * phi[n, j] - phi[n-1, j])
                    end
                elseif V[i, j] + V[i-1, j] < 0 && V[i, j] + V[i-1, j] < 0
                    if to_val[5] == true
                        val_push[5] += rho * V[i, j] * (-0.5)
                    else
                        b_push -= rho * V[i, j] * (-0.5) * (2 * phi[n, j] - phi[n-1, j])
                    end
                    val_push[6] += rho * V[i, j] * (2)
                    val_push[7] += rho * V[i, j] * (-1.5)
                elseif V[i, j] + V[i, j+1] <= 0 && V[i, j-1] + V[i, j] >= 0
                    val_push[6] += rho * V[i, j] * (+0.5)
                    val_push[8] += rho * V[i, j] * (-0.5)
                else
                    if to_val[5] == true
                        val_push[5] += rho * V[i, j] * (-0.5)
                    else
                        b_push -= rho * V[i, j] * (-0.5) * (2 * phi[n, j] - phi[n-1, j])
                    end
                    val_push[6] += rho * V[i, j] * (+1.5)
                    val_push[8] += rho * V[i, j] * (-1.5)
                    if to_val[9] == true
                        val_push[9] += rho * V[i, j] * (0.5)
                    else
                        b_push -= rho * V[i, j] * (0.5) * (2 * phi[n, j] - phi[n-1, j])
                    end
                end
                Ap[i,j]=val_push[7]
                val_push = val_push[findall(x -> x == true, to_val)]
                push!(val, val_push...)
                push!(b, b_push)
                continue
            end

            ## 处理最外层边界节点
            bd_coff = toABC(n,m,i, j,alpha,gamma,phi_symbol, bounds)
            if bd_coff[2] == 0
                val_push[7] = 1
                val_push = val_push[findall(x -> x == true, to_val)]
                push!(val, val_push...)
                push!(b, bd_coff[3] / bd_coff[1])
                continue
            end

            if i == 1 && j == 1
                val_push = 0.5*[
                    0,
                    -0.5 * rho * V[1, 1] - mu * 0.5 * (gamma[1, 1] / Ja[1, 1] + gamma[2, 1] / Ja[2, 1] + beta[1, 1] / Ja[1, 1] + beta[1, 2] / Ja[1, 2]),
                    0,
                    0.5 * rho * U[1, 1] - mu * 0.5 * (beta[1, 1] / Ja[1, 1] + beta[2, 1] / Ja[2, 1] + alpha[1, 1] / Ja[1, 1] + alpha[1, 2] / Ja[1, 2]),
                    0,
                    0
                ]
                val_push[1] = -val_push[2] - val_push[4] + bd_coff[1] / bd_coff[2]

                val_push[1]+=rho*v[i,j]/y_uv[i,j]*Ja[i,j]*0.25
                val_push[1]-=(mu*x_v[i,j]/y_uv[i,j]/2+mu*x_u[i,j]/y_uv[i,j]/2)*0.5
                val_push[2]+=mu*x_u[i,j]/y_uv[i,j]/2*0.5
                val_push[4]+=mu*x_v[i,j]/y_uv[i,j]/2*0.5

                push!(val, val_push...)
                push!(b, bd_coff[3] / bd_coff[2] + S[i, j] * 0.25)
                continue
            end
            if i == n && j == 1
                val_push = 0.5*[
                    0,
                    0.5 * rho * V[n, 1] - mu * 0.5 * (gamma[n, 1] / Ja[n, 1] + gamma[n-1, 1] / Ja[n-1, 1] + beta[n, 1] / Ja[n, 1] + beta[n, 2] / Ja[n, 2]),
                    0,
                    0,
                    0.5 * rho * U[n, 1] - mu * 0.5 * (beta[n, 1] / Ja[n, 1] + beta[n-1, 1] / Ja[n-1, 1] + alpha[n, 1] / Ja[n, 1] + alpha[n, 2] / Ja[n, 2]),
                    0
                ]
                val_push[3] = -val_push[2] - val_push[5] + bd_coff[1] / bd_coff[2]

                push!(val, val_push...)
                push!(b, bd_coff[3] / bd_coff[2] + S[i, j] * 0.25)
                continue
            end
            if i == n && j == m
                val_push = 0.5*[
                    0,
                    0,
                    -0.5 * rho * U[n, m] - mu * 0.5 * (beta[n, m] / Ja[n, m] + beta[n-1, m] / Ja[n-1, m] + alpha[n, m-1] / Ja[n, m-1] + alpha[n, m] / Ja[n, m]),
                    0,
                    0.5 * rho * V[n, m] - mu * 0.5 * (gamma[n, m] / Ja[n, m] + gamma[n-1, m] / Ja[n-1, m] + beta[n, m] / Ja[n, m] + beta[n, m-1] / Ja[n, m-1]),
                    0
                ]
                val_push[6] = -val_push[3] - val_push[5] + bd_coff[1] / bd_coff[2]

                push!(val, val_push...)
                push!(b, bd_coff[3] / bd_coff[2] + S[i, j] * 0.25)
                continue
            end
            if i == 1 && j == m
                val_push = 0.5*[
                    0,
                    -0.5 * rho * U[1, m] - mu * 0.5 * (beta[1, m] / Ja[1, m] + beta[2, m] / Ja[2, m] + alpha[1, m-1] / Ja[1, m-1] + alpha[1, m] / Ja[1, m]),
                    0,
                    0,
                    0.5 * rho * V[n, m] - mu * 0.5 * (gamma[1, m] / Ja[1, m] + gamma[2, m] / Ja[2, m] + beta[1, m] / Ja[1, m] + beta[1, m-1] / Ja[1, m-1]),
                    0
                ]
                val_push[4] = -val_push[2] - val_push[5] + bd_coff[1] / bd_coff[2]

                val_push[4]+=rho*v[i,j]/y_uv[i,j]*Ja[i,j]*0.25
                val_push[4]+=(mu*x_v[i,j]/y_uv[i,j]/2-mu*x_u[i,j]/y_uv[i,j]/2)*0.5
                val_push[2]-=mu*x_v[i,j]/y_uv[i,j]/2*0.5
                val_push[5]+=mu*x_u[i,j]/y_uv[i,j]/2*0.5

                push!(val, val_push...)
                push!(b, bd_coff[3] / bd_coff[2] + S[i, j] * 0.25)
                continue
            end
            if i == 1
                val_push[8] = -rho * V[i, j] / 2 - mu * 0.25 * (beta[i, j+1] / Ja[i, j+1] - beta[i, j-1] / Ja[i, j-1]) - mu * 0.5 * (gamma[i, j] / Ja[i, j] + gamma[i-1, j] / Ja[i-1, j])
                val_push[11] = +rho * U[i, j] / 4 - mu * 0.25 * (alpha[i, j] / Ja[i, j] + alpha[i, j+1] / Ja[i, j+1]) - mu * 0.25 * (beta[i, j] / Ja[i, j] + beta[i+1, j] / Ja[i+1, j])
                val_push[3] = -rho * U[i, j] / 4 - mu * 0.25 * (alpha[i, j] / Ja[i, j] + alpha[i, j-1] / Ja[i, j-1]) + mu * 0.25 * (beta[i, j] / Ja[i, j] + beta[i-1, j] / Ja[i-1, j])
                val_push[7] = -sum(val_push) + bd_coff[1] / bd_coff[2] * sqrt(gamma[i, j])
                

                val_push[7]+=rho*v[i,j]/y_uv[i,j]*Ja[i,j]*0.5
                val_push[3]-=mu*x_v[i,j]/y_uv[i,j]/2*0.5
                val_push[7]-=mu*x_u[i,j]/y_uv[i,j]/2
                val_push[8]+=mu*x_u[i,j]/y_uv[i,j]/2
                val_push[11]+=mu*x_v[i,j]/y_uv[i,j]/2*0.5


                push!(val, val_push[findall(x -> x == true, to_val)]...)
                push!(b, bd_coff[3] / bd_coff[2] * sqrt(gamma[i, j]) + S[i, j] * 0.5)
                continue
            end
            if i == n
                val_push[6] = rho * V[i, j] / 2 + mu * 0.25 * (beta[i, j+1] / Ja[i, j+1] - beta[i, j-1] / Ja[i, j-1]) - mu * 0.5 * (gamma[i, j] / Ja[i, j] + gamma[i-1, j] / Ja[i-1, j])
                val_push[11] = +rho * U[i, j] / 4 - mu * 0.25 * (alpha[i, j] / Ja[i, j] + alpha[i, j+1] / Ja[i, j+1]) + mu * 0.25 * (beta[i, j] / Ja[i, j] + beta[i-1, j] / Ja[i-1, j])
                val_push[3] = -rho * U[i, j] / 4 - mu * 0.25 * (alpha[i, j] / Ja[i, j] + alpha[i, j-1] / Ja[i, j-1]) - mu * 0.25 * (beta[i, j] / Ja[i, j] + beta[i-1, j] / Ja[i-1, j])
                val_push[7] = -sum(val_push) + bd_coff[1] / bd_coff[2] * sqrt(gamma[i, j])

                push!(val, val_push[findall(x -> x == true, to_val)]...)
                push!(b, bd_coff[3] / bd_coff[2] * sqrt(gamma[i, j]) + S[i, j] * 0.5)
                continue
            end
            if j == 1
                val_push[6] = rho * V[i, j] / 4 + mu * 0.25 * (beta[i, j] / Ja[i, j] + beta[i, j+1] / Ja[i, j+1]) - mu * 0.25 * (gamma[i, j] / Ja[i, j] + gamma[i-1, j] / Ja[i-1, j])
                val_push[8] = -rho * V[i, j] / 4 - mu * 0.25 * (beta[i, j] / Ja[i, j] + beta[i, j+1] / Ja[i, j+1]) - mu * 0.25 * (gamma[i, j] / Ja[i, j] + gamma[i+1, j] / Ja[i+1, j])
                val_push[11] = rho * U[i, j] / 2 - mu * 0.5 * (alpha[i, j] / Ja[i, j] + alpha[i, j+1] / Ja[i, j+1]) + mu * 0.25 * (beta[i-1, j] / Ja[i-1, j] - beta[i+1, j] / Ja[i+1, j])
                val_push[7] = -sum(val_push) + bd_coff[1] / bd_coff[2] * sqrt(alpha[i, j])


                val_push[7]+=rho*v[i,j]/y_uv[i,j]*Ja[i,j]*0.5
                val_push[7]-=mu*x_v[i,j]/y_uv[i,j]/2
                val_push[6]-=mu*x_u[i,j]/y_uv[i,j]/2*0.5
                val_push[8]+=mu*x_u[i,j]/y_uv[i,j]/2*0.5
                val_push[11]+=mu*x_v[i,j]/y_uv[i,j]/2

                push!(val, val_push[findall(x -> x == true, to_val)]...)
                push!(b, bd_coff[3] / bd_coff[2] * sqrt(alpha[i, j]) + S[i, j] * 0.5)
                continue
            end
            if j == m
                val_push[6] = rho * V[i, j] / 4 - mu * 0.25 * (beta[i, j] / Ja[i, j] + beta[i, j-1] / Ja[i, j-1]) - mu * 0.25 * (gamma[i, j] / Ja[i, j] + gamma[i-1, j] / Ja[i-1, j])
                val_push[8] = -rho * V[i, j] / 4 + mu * 0.25 * (beta[i, j] / Ja[i, j] + beta[i, j-1] / Ja[i, j-1]) - mu * 0.25 * (gamma[i, j] / Ja[i, j] + gamma[i+1, j] / Ja[i+1, j])
                val_push[3] = -rho * U[i, j] / 2 - mu * 0.5 * (alpha[i, j] / Ja[i, j] + alpha[i, j-1] / Ja[i, j-1]) - mu * 0.25 * (beta[i-1, j] / Ja[i-1, j] - beta[i+1, j] / Ja[i+1, j])
                val_push[7] = -sum(val_push) + bd_coff[1] / bd_coff[2] * sqrt(alpha[i, j])

                val_push[7]+=rho*v[i,j]/y_uv[i,j]*Ja[i,j]*0.5
                val_push[3]-=mu*x_v[i,j]/y_uv[i,j]/2
                val_push[6]-=mu*x_u[i,j]/y_uv[i,j]/2*0.5
                val_push[8]+=mu*x_u[i,j]/y_uv[i,j]/2*0.5
                val_push[7]+=mu*x_v[i,j]/y_uv[i,j]/2

                push!(val, val_push[findall(x -> x == true, to_val)]...)
                push!(b, bd_coff[3] / bd_coff[2] * sqrt(alpha[i, j]) + S[i, j] * 0.5)
                continue
            end
        end
    end
    Ap[2:end-1,1]=Ap[2:end-1,2]
    Ap[2:end-1,m]=Ap[2:end-1,m-1]
    Ap[1,2:end-1]=Ap[2,2:end-1]
    Ap[n,2:end-1]=Ap[n-1,2:end-1]
end
