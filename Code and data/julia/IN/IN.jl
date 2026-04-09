using Lux, DiffEqFlux, DifferentialEquations, Optimization, OptimizationOptimJL, Random, Plots
using DataFrames, CSV, ComponentArrays, OptimizationOptimisers
using FiniteDiff, LinearAlgebra, Distributions, Statistics

rng = Random.default_rng()
Random.seed!(14);
const SAVE_PATH = "C:/Users/XLYC/Desktop/流感代码修改/python/IN/"

# ============================================================
# 数据准备
# ============================================================
source_data = DataFrame(CSV.File("C:/Users/XLYC/Desktop/流感代码修改/julia/IN/total/INEMA.csv"))
n = 0
m = 312
data_weekly = source_data.EMA[(n+1):n+m]
display(plot(data_weekly, label="Weekly Cases", lw=2))
println(length(data_weekly))
data_acc = zeros(Float32, length(data_weekly))
data_acc[1] = data_weekly[1]
for j in 2:length(data_weekly)
    data_acc[j] = sum(data_weekly[1:j])
end
trainingdata_1 = Float32.(data_acc)
trainingdata_2 = Float32.(data_weekly)
N = 6630000

# ============================================================
# 神经网络与模型
# ============================================================
ann = Lux.Chain(
    Lux.Dense(3, 16, tanh),
    Lux.Dense(16, 1, softplus)
)
p_0, st = Lux.setup(rng, ann)

function SEIR_nn(du, u, p, t)
    S, E, I, R, C = u
    t1 = sin(t * 2.0f0 * Float32(π) / 52.0f0)
    t2 = -cos(t * 2.0f0 * Float32(π) / 52.0f0) + 1
    t3 = t/156
    beta = min((ann([ t1, t2,t3], p, st)[1])[1], 50.0)
    du[1] = -beta * S * I / N + 0.0078 * R 
    du[2] = beta * S * I / N - 3.5 * E
    du[3] = 3.5 * E - I
    du[4] = I - 0.0078 * R
    du[5] = 3.5 * E
end

u0 = Float32[N-57, 13, 44, 0, 44]
tspan = (0.0, 311.0)
tsteps = range(tspan[1], tspan[2], length = length(data_weekly))
prob_neuralode = ODEProblem(SEIR_nn, u0, tspan, ComponentArray(p_0))

days = 312 * 7
tsteps_daily = range(tspan[1], tspan[2], length = days)

# ============================================================
# 预测函数
# ============================================================
function predict_neuralode(θ)
    prob = remake(prob_neuralode, p = θ)
    sol = solve(prob, Tsit5(), saveat = tsteps, abstol = 1e-8, reltol = 1e-6)
    
    if sol.retcode != :Success
        return nothing
    end
    
    sol_arr = Array(sol)
    C = sol_arr[5, :]
    D = vcat(max(C[1], 1.0), max.(diff(C), 1.0))
    
    return C, D
end

# ============================================================
# 修改后的损失函数 - 添加第一年约束
# ============================================================
function loss_neuralode(p)
    result = predict_neuralode(p)
    if result === nothing
        return 1e10f0
    end
    
    pred_C, pred_D = result
    y = trainingdata_1
    loss = sum(abs2, log.(y) .- log.(pred_C))
    # λ = max.(pred_D, 1.0f0)
    
    # # 基础损失
    # poisson_nll = sum(λ .- y .* log.(λ))
    
    # # ========== Top-K 高峰损失 ==========
    # K = 5
    # sorted_indices = sortperm(y, rev=true)
    # top_k_indices = sorted_indices[1:K]
    # top_k_pred = pred_D[top_k_indices]
    # top_k_obs = y[top_k_indices]
    # top_k_loss = sum(((top_k_pred .- top_k_obs) ./ top_k_obs).^2)
    # max_loss = abs(maximum(pred_D) - maximum(y)) / maximum(y)
    
    # # ========== 第一年约束损失（新增）==========
    # first_year_end = 52  # 第一年52周
    # first_year_pred = pred_D[1:first_year_end]
    # first_year_obs = y[1:first_year_end]
    
    # # 方法1: 不对称惩罚 - 仅惩罚预测值高于实际值的情况
    # #overestimate_penalty = sum(max.(first_year_pred .- first_year_obs, 0.0f0).^2 ./ 
    #  #                          (first_year_obs.^2 .+ 1.0f0))
    
    # # 方法2: 第一年累计值约束
    # first_year_cum_pred = sum(first_year_pred)
    # first_year_cum_obs = sum(first_year_obs)
    # cum_loss = ((first_year_cum_pred - first_year_cum_obs) / first_year_cum_obs)^2
    
    # # 方法3: 第一年每周相对误差（带不对称权重）
    # #relative_errors = (first_year_pred .- first_year_obs) ./ (first_year_obs .+ 1.0f0)
    # # 对高估情况给更大权重
    # #asymmetric_weights = ifelse.(relative_errors .> 0, 2.0f0, 1.0f0)
    # #first_year_weighted_loss = sum(asymmetric_weights .* relative_errors.^2)
    
    #return poisson_nll + 2000.0f0 * top_k_loss + 1000.0f0 * max_loss + 2500.0f0 * cum_loss
           #3000.0f0 * overestimate_penalty + 2000.0f0 * cum_loss + 
           #1500.0f0 * first_year_weighted_loss
           return loss
end

# ============================================================
# 训练
# ============================================================
println("="^50)
println("开始训练...")
println("="^50)

loss_curve = Float64[]

callback = function (p, l; doplot = false)
    if isfinite(l)
        push!(loss_curve, l)
    end
    if length(loss_curve) % 50 == 0
        println("Iteration $(length(loss_curve)): Loss = $(round(l, digits=2))")
    end
    return false
end

pinit = ComponentArray(p_0)
adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss_neuralode(x), adtype)
optprob = Optimization.OptimizationProblem(optf, pinit)

# ADAM
result_1 = Optimization.solve(optprob, OptimizationOptimisers.ADAM(0.01),
    callback = callback, maxiters = 1000)

# LBFGS
optprob2 = remake(optprob, u0 = result_1.u)
pfinal = try
    result_2 = Optimization.solve(optprob2, Optim.LBFGS(),
        callback = callback, allow_f_increases = false)
    result_2.u
catch
    println("LBFGS 失败，使用 ADAM 结果")
    result_1.u
end

println("训练完成! 最终 Loss: ", round(loss_neuralode(pfinal), digits=2))

# ============================================================
# 计算每日 beta 值
# ============================================================
function compute_daily_beta(θ, days=312*7)
    beta_values = Float64[]
    for d in 0:(days-1)
        t = d / 7.0  # 转换为周
        t1 = sin(t * 2.0 * π / 52.0)
        t2 = -cos(t * 2.0 * π / 52.0) + 1
        t3 = t / 156.0
        beta = min((ann([Float32(t1), Float32(t2),Float32(t3)], θ, st)[1])[1], 50.0)
        push!(beta_values, beta)
    end
    return beta_values
end

# ============================================================
# 输出图像
# ============================================================
pred = predict_neuralode(pfinal)

# 图1: 每周新增病例对比图
plt1 = plot(1:312, trainingdata_2, 
    label="REAL", 
    lw=2, 
    color=:blue,
    marker=:circle,
    markersize=2,
    alpha=0.7)
plot!(plt1, 1:312, pred[2], 
    label="fit", 
    lw=2, 
    color=:red,
    linestyle=:dash)
# 标记第一年区域
xlabel!(plt1, "WEEK")
ylabel!(plt1, "New cases")
display(plt1)
savefig(plt1, SAVE_PATH * "weekly_cases.png")

# 图2: 每日 beta 值
daily_beta = compute_daily_beta(pfinal)
plt2 = plot(1:length(daily_beta), daily_beta,
    label=" β ",
    lw=1.5,
    color=:green,
    xlabel="day",
    ylabel="β ",
    title="(β)")
# 添加年份分隔线
for year in 1:6
    vline!(plt2, [year * 365], label="", color=:gray, linestyle=:dash, alpha=0.5)
end
display(plt2)
savefig(plt2, SAVE_PATH * "daily_beta.png")

println("图像已保存至: ", SAVE_PATH)
# ============================================================
# 拉普拉斯近似（简化版 - 仅协方差缩放）
# ============================================================
println("开始拉普拉斯近似...")

param_vec = Vector{Float64}(vec(pfinal))
loss_for_hessian(v) = Float64(loss_neuralode(ComponentArray(v, getaxes(pfinal))))

# 计算Hessian并正则化
H = FiniteDiff.finite_difference_hessian(loss_for_hessian, param_vec)
H = (H + H') / 2
min_eig = minimum(eigvals(H))
if min_eig <= 0
    H = H + (abs(min_eig) + 1.0) * I
end

# 计算协方差矩阵
Σ = inv(H)
Σ = (Σ + Σ') / 2

# 缩放协方差矩阵：限制相对标准差不超过20%
param_std = sqrt.(max.(diag(Σ), 0.0))
relative_std = param_std ./ (abs.(param_vec) .+ 1e-6)
max_allowed_std = 2.5
if maximum(relative_std) > max_allowed_std
    scale_factor = (max_allowed_std / maximum(relative_std))^2
    Σ = Σ * scale_factor
    println("协方差缩放因子: ", round(sqrt(scale_factor), digits=4))
end

# 确保正定
min_eig_cov = minimum(eigvals(Σ))
if min_eig_cov <= 0
    Σ = Σ + (abs(min_eig_cov) + 1e-8) * I
end

# 采样
n_samples = 500
posterior = MvNormal(param_vec, Σ)
param_samples = rand(posterior, n_samples)

# ============================================================
# 采样函数
# ============================================================
function sample_predictions(param_samples, pfinal)
    n_samples = size(param_samples, 2)
    weekly_preds = zeros(312, n_samples)
    valid = 0
    
    for i in 1:n_samples
        p = ComponentArray(param_samples[:, i], getaxes(pfinal))
        result = predict_neuralode(p)
        if result !== nothing
            _, D = result
            if all(isfinite.(D)) && all(D .> 0) && all(D .< 1e7)
                valid += 1
                weekly_preds[:, valid] = D
            end
        end
    end
    return weekly_preds[:, 1:valid], valid
end

function sample_beta(param_samples, pfinal, n_days)
    n_samples = size(param_samples, 2)
    beta_preds = zeros(n_days, n_samples)
    valid = 0
    
    for i in 1:n_samples
        p = ComponentArray(param_samples[:, i], getaxes(pfinal))
        betas = zeros(n_days)
        ok = true
        
        for d in 0:(n_days-1)
            t = d / 7.0
            t1 = Float32(sin(t * 2π / 52.0) )
            t2 = Float32(-cos(t * 2π / 52.0) + 1)
            t3 = Float32(t / 156.0)
            b = (ann([t1, t2, t3], p, st)[1])[1]
            if !isfinite(b) || b < 0 || b > 100
                ok = false
                break
            end
            betas[d+1] = min(b, 50.0)
        end
        
        if ok
            valid += 1
            beta_preds[:, valid] = betas
        end
    end
    return beta_preds[:, 1:valid], valid
end

function get_optimal_beta(pfinal, n_days)
    betas = zeros(n_days)
    for d in 0:(n_days-1)
        t = d / 7.0
        t1 = Float32(sin(t * 2π / 52.0))
        t2 = Float32(-cos(t * 2π / 52.0) + 1)
        t3 = Float32(t / 156.0)
        betas[d+1] = min((ann([t1, t2, t3], pfinal, st)[1])[1], 50.0)
    end
    return betas
end

# ============================================================
# 计算预测和置信区间
# ============================================================
println("计算置信区间...")
n_days = 312 * 7

weekly_preds, n_valid_weekly = sample_predictions(param_samples, pfinal)
beta_preds, n_valid_beta = sample_beta(param_samples, pfinal, n_days)

println("有效样本: weekly=$n_valid_weekly, beta=$n_valid_beta")

# 最优预测
pred_optimal = predict_neuralode(pfinal)
weekly_optimal = pred_optimal[2]
beta_optimal = get_optimal_beta(pfinal, n_days)

# 95%置信区间
weekly_lower = [quantile(weekly_preds[i, :], 0.025) for i in 1:312]
weekly_upper = [quantile(weekly_preds[i, :], 0.975) for i in 1:312]
beta_lower = [quantile(beta_preds[i, :], 0.025) for i in 1:n_days]
beta_upper = [quantile(beta_preds[i, :], 0.975) for i in 1:n_days]

# ============================================================
# 绘图
# ============================================================
# 图1: 每周新增病例
plt1 = plot(size=(1000, 600), dpi=150)
plot!(plt1, 1:312, weekly_lower, fillrange=weekly_upper, fillalpha=0.3, color=:red, label="95% CI", linealpha=0)
plot!(plt1, 1:312, weekly_optimal, lw=2, color=:red, label="fit")
scatter!(plt1, 1:312, trainingdata_2, color=:blue, markersize=3, alpha=0.7, label="real")
xlabel!(plt1, "week"); ylabel!(plt1, "New cases");
display(plt1)
savefig(plt1, SAVE_PATH * "weekly_cases_with_CI.png")

# 图2: 每日beta
plt2 = plot(size=(1000, 600), dpi=150)
plot!(plt2, 1:n_days, beta_lower, fillrange=beta_upper, fillalpha=0.3, color=:green, label="95% CI", linealpha=0)
plot!(plt2, 1:n_days, beta_optimal, lw=1.5, color=:darkgreen, label="β ")
for year in 1:6
    vline!(plt2, [year * 365], label="", color=:gray, linestyle=:dash, alpha=0.5)
end
xlabel!(plt2, "Day"); ylabel!(plt2, "β "); title!(plt2, " β")
display(plt2)
savefig(plt2, SAVE_PATH * "daily_beta_with_CI.png")

# ============================================================
# 保存CSV
# ============================================================
weekly_df = DataFrame(
    Week = 1:312,
    Observed = trainingdata_2,
    Predicted = weekly_optimal,
    CI_Lower = weekly_lower,
    CI_Upper = weekly_upper
)
CSV.write(SAVE_PATH * "weekly_cases_results.csv", weekly_df)

beta_df = DataFrame(
    Day = 1:n_days,
    Beta = beta_optimal,
    CI_Lower = beta_lower,
    CI_Upper = beta_upper
)
CSV.write(SAVE_PATH * "daily_beta_results.csv", beta_df)

println("完成! 文件已保存至: ", SAVE_PATH)