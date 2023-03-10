module EasyCFD2D

using SparseArrays
#using EquationsSolver

include("./BaseFunctions.jl")
#include("./PostProcess.jl")
include("./Grider/CofficientGeneration.jl")
include("./Grider/GriderMethods.jl")

#边界数据结构
export bound

#网格划分方法
export GSGrider, JacobianGrider

#作图方法
#export gridPlot

end # module EasyCFD2D
