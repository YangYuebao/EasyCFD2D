abstract type HFProblem end

struct FluidProperties
    rho::Float64
    mu::Float64
end

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

struct FluidProblem <: HFProblem
    bounds::Vector{bound}
    grid::Grid
    property::FluidProperties
    coodinate::AbstractCoodinateTypes
end