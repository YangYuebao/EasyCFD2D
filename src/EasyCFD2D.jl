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
include("./FlowSolver/Solver.jl")

#边界数据结构
export bound

#网格划分方法
export GSGrider, JacobianGrider

# 边界类型
export stillWall,movingWall,FDOutlet,
    velocityInlet,pressureInlet,symetryAxis,
    generalBound

# 坐标系类型
export Rectangular,Cylindrical

# 对流离散格式
export SecondOrderUpwind

# 求解流场
export solvefield

# 压力速度耦合算法
export SIMPLE

# test
export to_val_index,toABC,getUV,fieldDiff

end # module EasyCFD2D
