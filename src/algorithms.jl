immutable FEMDiffEqPoisson <: AbstractPoissonFEMAlgorithm end

immutable FEMDiffEqHeatEuler <: AbstractHeatFEMAlgorithm end
immutable FEMDiffEqHeatImplicitEuler <: AbstractHeatFEMAlgorithm end
immutable FEMDiffEqHeatCrankNicholson <: AbstractHeatFEMAlgorithm end
immutable FEMDiffEqHeatSemiImplicitEuler <: AbstractHeatFEMAlgorithm end
immutable FEMDiffEqHeatSemiImplicitCrankNicholson <: AbstractHeatFEMAlgorithm end
