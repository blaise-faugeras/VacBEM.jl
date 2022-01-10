module Assembly

using LinearAlgebra
using SparseArrays

using GreenDeltaStar

include("assemPotentialKandV.jl")
export assemPotentialKandV

include("assemPotentialKandVwithDuffy.jl")
export assemPotentialKandVwithDuffy

include("evalPsiFromBEM.jl")
export evalPsiFromBEM

include("rectangleCoil.jl")
export evalPsiFromRectangleCoil
export assemblyRhsRectangleCoil

end
