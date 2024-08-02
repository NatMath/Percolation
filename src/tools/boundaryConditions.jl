abstract type BoundaryCondition end
struct NormalBoundary <: BoundaryCondition end
export NormalBoundary
struct PeriodicBoundary <: BoundaryCondition end
export PeriodicBoundary