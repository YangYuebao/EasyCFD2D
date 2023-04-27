abstract type HFProblem end

# 流体参数
struct FluidProperties
    rho::Float64
    mu::Float64
end

# 网格信息
struct Grid
    x_uv::Matrix{Float64}
    y_uv::Matrix{Float64}
    x_u::Matrix{Float64}
    x_v::Matrix{Float64}
    y_u::Matrix{Float64}
    y_v::Matrix{Float64}
    alpha::Matrix{Float64}
    beta::Matrix{Float64}
    gamma::Matrix{Float64}
    Ja::Matrix{Float64}
end

# 问题定义
struct FluidProblem <: HFProblem
    bounds::Vector{bound}
    grid::Grid
    property::FluidProperties
    coodinate::AbstractCoodinateTypes
end