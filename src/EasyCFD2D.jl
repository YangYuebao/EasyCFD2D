module EasyCFD2D

using SparseArrays
#using EquationsSolver

include("./BasicTypes/BoundaryTypes.jl")
include("./BasicTypes/CoodinateTypes.jl")
include("./BasicTypes/ProblemTypes.jl")
include("./BasicTypes/BaseFunctions.jl")

include("./Grider/CofficientGeneration.jl")
include("./Grider/GriderMethods.jl")

include("./FlowSolver/FieldInitialize.jl")
include("./FlowSolver/FlowTerm.jl")
include("./FlowSolver/PressureVelocityCouple.jl")
include("./FlowSolver/RenewCoff.jl")
include("./FlowSolver/SourceTerm.jl")

#边界数据结构
export bound

#网格划分方法
export GSGrider, JacobianGrider

export stillWall,movingWall,FDOutlet,
    velocityInlet,pressureInlet,symetryAxis,
    generalBound

# 对流离散格式
export SecondOrderUpwind

# test
export to_val_index,toABC,getUV,SIMPLE

end # module EasyCFD2D
