using Lux, DiffEqFlux, DifferentialEquations, Optimization, OptimizationOptimJL, Random, Plots
using DataFrames
using CSV
using ComponentArrays
using OptimizationOptimisers

rng = Random.default_rng()
Random.seed!(14);

source_data = DataFrame(CSV.File("C:/Users/Administrator/Desktop/Newproject/data/iowa2016-2023.csv"))
n = 0
m = 377
data_weekly = source_data.ILITotal[(n+1):n+m+1]
display(plot(data_weekly, label="Weekly Cases", lw=2))
println(length(data_weekly))
trainingdata = Float32.(data_weekly)

ann = Lux.Chain(Lux.Dense(1, 32, tanh), Lux.Dense(32, 1))
p_0, st = Lux.setup(rng, ann)
function SIR_nn(du, u, p, t)
    I, R = u
    
    du[1] = abs(((ann([t], p, st))[1])[1]) - I
    du[2] = I 
end
u0 = Float32[data_weekly[1], 0]
tspan = (0.0f0, 377.0f0)
tsteps = range(tspan[1], tspan[2], length=length(data_weekly))
prob_neuralode = ODEProblem(SIR_nn, u0, tspan, ComponentArray(p_0))



function predict_neuralode(θ)
    #Array(prob_neuralode(u0, p, st)[1])
    prob = remake(prob_neuralode, p=θ)
    Array(solve(prob, Tsit5(), saveat=tsteps))
    
    #Array(concrete_solve(prob, Tsit5(), u_0, θ, saveat=1,
        #abstol=1e-6, reltol=1e-6))
end




predict_neuralode(p_0)[1, :]
size(predict_neuralode(p_0)[1, :]) == size(data_weekly)

function loss_neuralode(p)
    pred = predict_neuralode(p)[1, :]
    loss =sum(abs2, log.(trainingdata) .- log.(pred))
    return loss, pred
end

loss_neuralode(p_0)

loss_curve=Float64[]
callback = function (p, l, pred; doplot=false)
    println(l)
    push!(loss_curve,l)
    if doplot
        plt = scatter(tsteps, trainingdata, label="Weekly cases")
        plot!(plt, tsteps, pred, label="Predicted cases")
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
    maxiters=415)


optprob2 = remake(optprob, u0=result_neuralode.u)

result_neuralode2 = Optimization.solve(optprob2,
    Optim.LBFGS(),
    callback=callback,
    allow_f_increases=false)


pfinal = result_neuralode2.u

callback(pfinal, loss_neuralode(pfinal)...; doplot=true)
display(plot(1:length(loss_curve), loss_curve, xlabel="Iteration", ylabel="Loss", label="Loss Curve"))
#plot(loss_neuralode(p))


using BSON: @save
@save "C:/Users/Administrator/Desktop/Newproject/saving/ann_nn_ir.bson" ann
psave = collect(pfinal)
@save "C:/Users/Administrator/Desktop/Newproject/saving/ann_para_irlbfgs.bson" psave
pred = predict_neuralode(pfinal)[1, :]
plt = scatter(tsteps, trainingdata, label="Weekly cases")
plot!(plt, tsteps, pred, label="Predicted cases")
display(plot(plt))
savefig("C:/Users/Administrator/Desktop/Newproject/saving/annepicase.png")


#for i in 1:120
    #t = tsteps[i]
    #a = ann([t],pfinal,st)[1][1]
    #b = abs(a)
    #c = pred[i]
    #println(b / c)
#end