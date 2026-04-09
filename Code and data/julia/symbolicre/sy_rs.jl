using SymbolicRegression
using Random: MersenneTwister
using Zygote
using MLJBase: machine, fit!, predict, report
using Test
using CSV
data = CSV.File("C:/Users/XLYC/Desktop/流感代码修改/R/SR/total_34nb.csv")

# 提取两列并转换为数值向量
col1 = data.temp  # 提取第一列
col2 = data.RR  # 提取第二列
x = Vector{Float64}(col1)
x= reshape(x, length(x), 1) 
y = Vector{Float64}(col2)

model = SRRegressor(
    binary_operators=[+, -, *, /],
    unary_operators=[exp],
    niterations=100
)
mach = machine(model, x, y)
fit!(mach)
