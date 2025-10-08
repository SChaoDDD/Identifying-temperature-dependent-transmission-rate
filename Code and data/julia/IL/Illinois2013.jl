using Lux, DiffEqFlux, DifferentialEquations, Optimization, OptimizationOptimJL, Random, Plots
using DataFrames
using CSV
using ComponentArrays
using OptimizationOptimisers
rng = Random.default_rng()
Random.seed!(14);

source_data = DataFrame(CSV.File("C:/Users/XLYX/Desktop/Illinois2013-2024_1.csv"))
n = 32
m = 54
data_weekly = source_data.EMA[(n+1):n+m+1]
display(plot(data_weekly, label="Weekly Cases", lw=2))
println(length(data_weekly))
data_acc = zeros(Float32, length(data_weekly))
data_acc[1] = data_weekly[1]
for j in 2:length(data_weekly)
    data_acc[j] = sum(data_weekly[1:j])  # 前 n 个元素的和
end
trainingdata_1 = Float32.(data_acc)
trainingdata_2 = Float32.(data_weekly)
N = 12900000



function simulate_f(t)
    return sin(0.1396*t)
end

ann = Lux.Chain(Lux.Dense(1, 32, tanh),  Lux.Dense(32, 1))
p_0, st = Lux.setup(rng, ann)
function SEIR_nn(du, u, p, t)
    S, E, I, R, C = u
    #total_population = S + E + I + R
    #scale_factor = N / total_population
    du[1] = - min(abs(((ann([t], p, st))[1])[1]), 5) * S * I / N + 0.0078 * R 
    du[2] = min(abs(((ann([t], p, st))[1])[1]), 5) * S * I / N - 3.5 * E
    du[3] = 3.5 * E - I
    du[4] = I - 0.0078 * R
    du[5] = 3.5 * E
    #du[1:4] .= du[1:4] .* scale_factor
end
u0 = Float32[12899538, 92, 370, 0, 370]
tspan = (0.0f0, 54.0f0)
tsteps = range(tspan[1], tspan[2], length=length(data_weekly))
prob_neuralode = ODEProblem(SEIR_nn, u0, tspan, ComponentArray(p_0))





function predict_neuralode(θ)
    #Array(prob_neuralode(u0, p, st)[1])
    prob = remake(prob_neuralode, p=θ)
    sol = Array(solve(prob, Tsit5(), saveat=tsteps))
    S = sol[1, :]
    E = sol[2, :]
    I = sol[3, :]
    R = sol[4, :]
    C = sol[5, :]
    D = zeros(Float32, length(C))
    C_1 = Float32.(C[1])
    C_2 = diff(C)
    D = vcat(C_1, C_2)
    return C, D
end




predict_neuralode(p_0)
size(predict_neuralode(p_0)[2]) == size(data_weekly)

function loss_neuralode(p)
    pred = predict_neuralode(p)
    loss = sum(abs2, log.(trainingdata_1) .- log.(pred[1]))
    return loss, pred
end

loss_neuralode(p_0)

loss_curve=Float64[]
callback = function (p, l, pred; doplot=false)
    println(l)
    push!(loss_curve,l)
    if doplot
        plt = scatter(tsteps, trainingdata_1, label="Accumulated Cases")
        plot!(plt, tsteps, pred[1], label="Predicted Accumulated Cases")
        display(plot(plt))
    end
    return false
end

pinit = ComponentArray(p_0)
callback(pinit, loss_neuralode(pinit)...; doplot=true)

adtype = Optimization.AutoZygote()

optf = Optimization.OptimizationFunction((x, p) -> loss_neuralode(x), adtype)
optprob = Optimization.OptimizationProblem(optf, pinit)

result_neuralode = Optimization.solve(optprob,
      OptimizationOptimisers.ADAM(0.01),
      callback=callback,
      maxiters=100)

optprob2 = remake(optprob, u0=result_neuralode.u)

result_neuralode2 = Optimization.solve(optprob2,
    Optim.LBFGS(),
    callback=callback,
    allow_f_increases=false)


pfinal = result_neuralode2.u

callback(pfinal, loss_neuralode(pfinal)...; doplot=true)
xticks = 1:10:length(loss_curve)
xticklabels = xticks
display(plot(1:length(loss_curve), loss_curve, xlabel="Iteration", ylabel="Loss", label="Loss Curve"))
  #plot(loss_neuralode(p))


  #using BSON: @save
  #@save "C:/Users/Administrator/Desktop/Newproject/saving/ann_nn_ir.bson_iowa" ann
  #psave = collect(pfinal)
  #@save "C:/Users/Administrator/Desktop/Newproject/saving/ann_para_irlbfgs.bson_iowa" psave
pred = predict_neuralode(pfinal)
plt = scatter(tsteps, trainingdata_1, label="Accumulated Cases")
plot!(plt, tsteps, pred[1], label="Predicted Accumulated Cases")
xlabel!(plt, "Week")
ylabel!(plt, "Cases")

display(plot(plt))

println("beta")
for i in 1:55
    t = tsteps[i]
    a = ann([t],pfinal,st)[1][1]
    b = abs(a)
    #c = pred[i]
    println(b)
end