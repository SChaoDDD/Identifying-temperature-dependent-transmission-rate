using SymbolicRegression
using Random: MersenneTwister
using Zygote
using MLJBase: machine, fit!, predict, report
using Test
using CSV
data = CSV.File("C:/Users/XLYC/Desktop/ILbeta.csv")


col1 = data.temp  
col2 = data.lag7  
x = Vector{Float64}(col1)
x= reshape(x, length(x), 1) 
y = Vector{Float64}(col2)

model = SRRegressor(
    binary_operators=[+, -, *, /],
    unary_operators=[exp],
    niterations=30
)
mach = machine(model, x, y)
fit!(mach)