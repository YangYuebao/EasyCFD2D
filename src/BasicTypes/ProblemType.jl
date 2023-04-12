abstract type HFProblem end

struct FluidProblem <: HFProblem
    bounds::Vector{bound}
    coodinate::AbstractCoodinateTypes
end