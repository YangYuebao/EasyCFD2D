# 源项暂时不处理
#=
function sorceTerm!(::Cylindrical, S::Matrix{Float64}, rho::Number, mu::Number, phi_field::Matrix{Float64}, v_field::Matrix{Float64}, y_uv::Matrix{Float64},x_u::Matrix{Float64},x_v::Matrix{Float64},Ja::Matrix{Float64})
    phi_u, phi_v = fieldDiff(phi_field)
    temp = ( - rho * v_field[1:end-1,:] .* phi_field[1:end-1,:] ./ y_uv[1:end-1,:]) .* Ja[1:end-1,:] + mu ./ y_uv[1:end-1,:] .* (-x_v[1:end-1,:] .* phi_u[1:end-1,:] + x_u[1:end-1,:] .* phi_v[1:end-1,:])
    S[1:end-1,:] += temp
    #S[end,:] += 2*temp[end,:] - temp[end-1,:]
    S
end

function sorceTerm!(::Rectangular, S::Matrix{Float64}, rho::Number, mu::Number, phi_field::Matrix{Float64}, v_field::Matrix{Float64}, y_uv::Matrix{Float64},x_u::Matrix{Float64},x_v::Matrix{Float64},Ja::Matrix{Float64})
    S
end
=#
