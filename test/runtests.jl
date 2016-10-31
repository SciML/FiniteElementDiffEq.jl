#!/usr/bin/env julia

const CPU_FLOPS = peakflops()
const TEST_PLOT = false
const LONGER_TESTS = true #Requires JLD
const TEST_CONDITIONAL_DEPS = true
const FILEIO_ENABLE = false
#Start Test Script
using FiniteElementDiffEq, DiffEqDevTools, DiffEqProblemLibrary
using Base.Test

tic()

# Internals
println("Quadrature Points Tests")
@time @test include("internals/quadpts_test.jl")
println("Assembly Tests")
@time @test include("internals/assembly_tests.jl")
println("Boundary Tests")
FILEIO_ENABLE && @time @test include("meshes/boundary_tests.jl")
println("Example Mesh Tests")
FILEIO_ENABLE && @time @test include("meshes/mesh_examples_tests.jl")
println("Simple Mesh Tests")
@time @test include("meshes/mesh_SimpleMesh_tests.jl")


#Finite Element
#Heat
println("Finite Element Heat Dt Tests")
@time @test include("heat/femheat_dtconvergence_tests.jl")
println("Finite Element Heat Dx Tests")
@time @test include("heat/femheat_dxconvergence_tests.jl")
println("Finite Element Heat Method Tests")
@time @test include("heat/femheat_methods_tests.jl")
println("Finite Element Nonlinear Heat Methods Tests")
@time @test include("heat/femheat_nonlinearmethods_tests.jl")
println("Finite Element Nonlinear System Heat Tests")
@time @test include("heat/femheat_system_tests.jl")
println("Heat Animation Test")
@time @test include("heat/femheat_animation_tests.jl")
println("Stochastic Heat Animation Test")
@time @test include("heat/femheat_stochasticanimation_tests.jl")

#Poisson
println("Finite Element Poisson Convergence Test")
@time @test include("poisson/fempoisson_convergence_tests.jl")
println("Finite Element Nonlinear Poisson Tests")
@time @test include("poisson/fempoisson_nonlinear_tests.jl")
println("Finite Element Nonlinear System Poisson Tests")
@time @test include("poisson/fempoisson_nonlinearsystem_tests.jl")
println("Finite Element Stochastic Poisson")
@time @test include("poisson/fempoisson_stochastic_tests.jl")
println("Finite Element Poisson")
@time @test include("poisson/fempoisson_linear_tests.jl")

println("Other Premades")
@time @test include("internals/other_premades_tests.jl")
println("Units Tests")
@time @test include("internals/units_tests.jl")

toc()
