using Pkg
Pkg.activate()
using SpecialFunctions
using Lux, DiffEqFlux, DifferentialEquations, Optimization, OptimizationOptimJL, Random, Plots
using DataFrames, CSV, ComponentArrays, OptimizationOptimisers
using LinearAlgebra, Distributions, Statistics
using Zygote, SciMLSensitivity

rng = Random.default_rng()
Random.seed!(14)
const SAVE_PATH = "C:/Users/XLYC/Desktop/流感代码修改/PLOS代码/JUlia/IA/IAuncertain1/"

# ============================================================
# 1. 数据准备
# ============================================================
source_data = DataFrame(CSV.File(
    "C:/Users/XLYC/Desktop/流感代码修改/PLOS代码/matlab/IAEMA.csv"))
n = 0
m = 312
data_weekly = source_data[(n+1):n+m, 1] * 2.74
display(plot(data_weekly, label="Weekly Cases", lw=2))
println("数据长度: ", length(data_weekly))

data_acc = zeros(Float32, length(data_weekly))
data_acc[1] = data_weekly[1]
for j in 2:length(data_weekly)
    data_acc[j] = sum(data_weekly[1:j])
end

trainingdata_1 = Float32.(data_acc)
trainingdata_2 = Float32.(data_weekly)
N = 3125800.0f0

# ============================================================
# 2. 网络构建函数
# ============================================================
function build_ann(seed::Int)
    local_rng = Random.MersenneTwister(seed)
    net = Lux.Chain(
        Lux.Dense(3, 16, tanh),
        Lux.Dense(16, 16, tanh),
        Lux.Dense(16, 1, softplus)
    )
    p, st = Lux.setup(local_rng, net)
    return net, p, st
end

# ============================================================
# 3. ODE 构建
# ============================================================
function make_SEIR_nn(net, st_fixed)
    function SEIR_nn(du, u, p, t)
        S, E, I, R, C = u
        t1 = sin(t * 2.0f0 * Float32(π) / 52.0f0)
        t2 = -cos(t * 2.0f0 * Float32(π) / 52.0f0) + 1.0f0
        t3 = t / 156.0f0
        nn_input = [t1, t2, t3]
        beta_raw = net(nn_input, p, st_fixed)[1][1]
        beta = clamp(beta_raw, 0.01f0, 50.0f0)
        S_ = max(S, 0.0f0)
        E_ = max(E, 0.0f0)
        I_ = max(I, 0.0f0)
        R_ = max(R, 0.0f0)
        du[1] = -beta * S_ * I_ / N + 0.0078f0 * R_
        du[2] =  beta * S_ * I_ / N - 3.5f0 * E_
        du[3] =  3.5f0 * E_ - I_
        du[4] =  I_ - 0.0078f0 * R_
        du[5] =  3.5f0 * E_
    end
    return SEIR_nn
end

u0    = Float32[N - 9, 2, 7, 0, 7]
tspan = (0.0f0, 311.0f0)
tsteps = range(tspan[1], tspan[2], length=312)
const SENSE_ALG = InterpolatingAdjoint(autojacvec=ZygoteVJP())

# ============================================================
# 4. 预测 / 损失 / Beta 函数生成器
# ============================================================
function make_predict(net, st_fixed)
    seir_fn   = make_SEIR_nn(net, st_fixed)
    p_dummy,_ = Lux.setup(Random.default_rng(), net)
    prob_base = ODEProblem(seir_fn, u0, tspan, ComponentArray(p_dummy))

    function predict_fn(θ)
        prob = remake(prob_base, p=θ)
        sol = solve(prob, Tsit5();
            saveat=tsteps, abstol=1e-6, reltol=1e-4,
            maxiters=100_000, verbose=false, sensealg=SENSE_ALG)
        if !SciMLBase.successful_retcode(sol) || size(Array(sol),2) != length(tsteps)
            return nothing
        end
        sol_arr = Array(sol)
        C = max.(sol_arr[5, :], 1.0f0)
        D = vcat(C[1], max.(diff(C), 1.0f0))
        return C, D
    end
    return predict_fn
end

function negbin_loglik(y::Real, μ::Real, r::Real)
    μ_d = max(Float64(μ), 1e-8)
    r_d = max(Float64(r), 1e-8)
    y_d = max(Float64(y), 0.0)
    ll  = loggamma(y_d + r_d) - loggamma(r_d) - loggamma(y_d + 1.0) +
          r_d * log(r_d / (r_d + μ_d)) + y_d * log(μ_d / (r_d + μ_d))
    return Float32(isfinite(ll) ? ll : -1f10)
end

function make_loss(predict_fn)
    function loss_negbin(p_full)
        nn_params = p_full.nn_params
        r = exp(Float32(p_full.log_r))
        result = predict_fn(nn_params)
        result === nothing && return 1.0f10
        _, pred_D = result
        nll = 0.0f0
        @inbounds for i in eachindex(trainingdata_2)
            nll -= negbin_loglik(trainingdata_2[i], pred_D[i], r)
        end
        return isfinite(nll) ? nll : 1.0f10
    end
    return loss_negbin
end

function make_beta_fn(net, st_fixed)
    function compute_daily_beta(θ, days=312*7)
        betas = Float64[]
        for d in 0:(days-1)
            t  = d / 7.0
            t1 = Float32(sin(t * 2π / 52))
            t2 = Float32(-cos(t * 2π / 52) + 1)
            t3 = Float32(t / 156)
            β  = clamp(Float64(net([t1, t2, t3], θ, st_fixed)[1][1]), 0.01, 50.0)
            push!(betas, β)
        end
        return betas
    end
    return compute_daily_beta
end

# ============================================================
# 5. 单成员训练
# ============================================================
function train_single_member(seed::Int, id::Int)
    println("\n" * "="^60)
    println("🔄 训练集成成员 #$id  (seed=$seed)")
    println("="^60)

    local_net, local_p0, local_st = build_ann(seed)
    predict_fn = make_predict(local_net, local_st)
    loss_fn    = make_loss(predict_fn)
    beta_fn    = make_beta_fn(local_net, local_st)

    p_init = ComponentArray(nn_params=ComponentArray(local_p0),
                            log_r=Float32(log(10.0)))
    println("  初始 NLL = ", loss_fn(p_init))

    cb = let cnt = Ref(0)
        function (state, l)
            isfinite(l) && (cnt[] += 1)
            cnt[] % 500 == 0 && println("  [#$id] Iter $(cnt[]): NLL=$(round(l,digits=2))")
            return false
        end
    end

    adtype  = Optimization.AutoZygote()
    optf    = Optimization.OptimizationFunction((x,_)->loss_fn(x), adtype)
    optprob = Optimization.OptimizationProblem(optf, p_init)

    result_adam = Optimization.solve(optprob,
        OptimizationOptimisers.ADAM(0.001), callback=cb, maxiters=3000)
    println("  [#$id] ADAM  NLL=$(round(loss_fn(result_adam.u),digits=2))")

    pfinal = try
        optprob2 = remake(optprob, u0=result_adam.u)
        res = Optimization.solve(optprob2, Optim.LBFGS();
            callback=cb, allow_f_increases=false, maxiters=3000)
        println("  [#$id] LBFGS NLL=$(round(loss_fn(res.u),digits=2))")
        res.u
    catch e
        println("  [#$id] LBFGS 失败 → 用 ADAM 结果")
        result_adam.u
    end

    final_nll = loss_fn(pfinal)
    r_val = exp(pfinal.log_r)
    println("  [#$id] ✅ 完成  NLL=$(round(final_nll,digits=2))  r=$(round(r_val,digits=3))")

    return (params=pfinal, predict_fn=predict_fn, beta_fn=beta_fn,
            nll=final_nll, r=r_val)
end

# ============================================================
# 6. 深度集成：训练全部成员
# ============================================================
const N_ENSEMBLE    = 15
const ENSEMBLE_SEEDS = [6,14,38, 42, 77,666, 123,347, 256, 512, 999,456,88,66,100]

println("\n🚀 开始深度集成训练 ($N_ENSEMBLE 个成员)")
members = [train_single_member(ENSEMBLE_SEEDS[i], i) for i in 1:N_ENSEMBLE]

# ============================================================
# 7. 筛选：选最优基准 + 剔除差模型
# ============================================================
nll_values = [m.nll for m in members]
println("\n" * "="^60)
println("📋 各成员 NLL 汇总")
println("-"^40)
for i in 1:N_ENSEMBLE
    println("  成员#$i  NLL=$(round(nll_values[i],digits=2))  r=$(round(members[i].r,digits=3))")
end

# —— 选出最佳模型作为基准 ——
best_idx = argmin(nll_values)
best_nll = nll_values[best_idx]
println("\n⭐ 最佳模型: 成员#$best_idx  NLL=$(round(best_nll,digits=2))")

# —— 剔除规则：NLL 超过最佳值 + 阈值 的视为明显差 ——
# 阈值取所有有效 NLL（排除 1e10）的标准差，或固定百分比
finite_nlls = filter(x -> x < 1.0f9, nll_values)
if length(finite_nlls) >= 3
    nll_median = median(finite_nlls)
    nll_iqr    = quantile(finite_nlls, 0.75) - quantile(finite_nlls, 0.25)
    # 保守规则：NLL > median + 1.5*IQR 则剔除（类似箱线图异常值）
    nll_cutoff = nll_median + 1.5 * max(nll_iqr, 0.01 * abs(nll_median))
else
    # 成员太少，放宽为最佳值的 10%
    nll_cutoff = best_nll * 1.10
end
println("  剔除阈值 NLL > $(round(nll_cutoff,digits=2))")

keep_mask = [(m.nll < 1.0f9) && (m.nll <= nll_cutoff)&&(m.r>3) for m in members]
for i in 1:N_ENSEMBLE
    tag = keep_mask[i] ? "✅ 保留" : "❌ 剔除"
    extra = (i == best_idx) ? " ⭐基准" : ""
    println("  成员#$i  NLL=$(round(nll_values[i],digits=2))  $tag$extra")
end

kept_indices = findall(keep_mask)
n_kept = length(kept_indices)
println("\n保留 $n_kept / $N_ENSEMBLE 个成员用于 Ensemble 不确定性估计")

if n_kept < 2
    println("⚠️ 保留成员不足 2 个，放宽阈值，保留所有有限 NLL 成员")
    keep_mask = [m.nll < 1.0f9 for m in members]
    kept_indices = findall(keep_mask)
    n_kept = length(kept_indices)
    println("  放宽后保留 $n_kept 个成员")
end

# ============================================================
# 8. 集成预测 + 不确定性
# ============================================================
println("\n计算集成预测...")
n_weeks = 312
n_days  = 312 * 7

# —— 新增病例 ——
all_pred_D = Matrix{Float64}(undef, n_weeks, n_kept)
all_pred_C = Matrix{Float64}(undef, n_weeks, n_kept)

for (col, idx) in enumerate(kept_indices)
    res = members[idx].predict_fn(members[idx].params.nn_params)
    if res !== nothing
        C_i, D_i = res
        all_pred_C[:, col] = Float64.(C_i)
        all_pred_D[:, col] = Float64.(D_i)
    else
        all_pred_C[:, col] .= NaN
        all_pred_D[:, col] .= NaN
        println("⚠️ 保留成员#$idx 预测失败")
    end
end

# 再次过滤 NaN
valid_cols_D = [!any(isnan.(all_pred_D[:, c])) for c in 1:n_kept]
all_pred_D = all_pred_D[:, valid_cols_D]
all_pred_C = all_pred_C[:, valid_cols_D]
n_final = size(all_pred_D, 2)
println("最终有效 Ensemble 成员数: $n_final")

ensemble_D_mean  = vec(mean(all_pred_D, dims=2))
ensemble_D_std   = vec(std(all_pred_D, dims=2))
ensemble_D_lower = ensemble_D_mean .- 1.96 .* ensemble_D_std
ensemble_D_upper = ensemble_D_mean .+ 1.96 .* ensemble_D_std
ensemble_D_q025  = vec(mapslices(x -> quantile(x, 0.025), all_pred_D, dims=2))
ensemble_D_q975  = vec(mapslices(x -> quantile(x, 0.975), all_pred_D, dims=2))

# —— Beta ——
all_betas = Matrix{Float64}(undef, n_days, n_final)
final_kept = kept_indices[valid_cols_D]   # 映射回原始成员索引
for (col, idx) in enumerate(final_kept)
    all_betas[:, col] = members[idx].beta_fn(members[idx].params.nn_params, n_days)
end

ensemble_beta_mean  = vec(mean(all_betas, dims=2))
ensemble_beta_std   = vec(std(all_betas, dims=2))
ensemble_beta_lower = max.(ensemble_beta_mean .- 1.96 .* ensemble_beta_std, 0.01)
ensemble_beta_upper = min.(ensemble_beta_mean .+ 1.96 .* ensemble_beta_std, 50.0)
ensemble_beta_q025  = max.(vec(mapslices(x -> quantile(x, 0.025), all_betas, dims=2)), 0.01)
ensemble_beta_q975  = min.(vec(mapslices(x -> quantile(x, 0.975), all_betas, dims=2)), 50.0)

# —— 最佳基准模型的预测（单独提取用于对比）——
best_member = members[best_idx]
best_result = best_member.predict_fn(best_member.params.nn_params)
best_C, best_D = best_result
best_beta = best_member.beta_fn(best_member.params.nn_params, n_days)

# ============================================================
# 9. 保存 CSV
# ============================================================
println("\n保存 CSV...")

df_beta = DataFrame(
    Day            = 1:n_days,
    Beta_Best      = best_beta,
    Beta_EnsMean   = ensemble_beta_mean,
    Beta_EnsStd    = ensemble_beta_std,
    Beta_Lower95   = ensemble_beta_lower,
    Beta_Upper95   = ensemble_beta_upper,
    Beta_Q025      = ensemble_beta_q025,
    Beta_Q975      = ensemble_beta_q975
)
CSV.write(SAVE_PATH * "daily_beta_ensemble.csv", df_beta)
println("✅ daily_beta_ensemble.csv  ($(nrow(df_beta)) 行)")

df_weekly = DataFrame(
    Week           = 1:n_weeks,
    Observed       = Float64.(trainingdata_2),
    Pred_Best      = Float64.(best_D),
    Pred_EnsMean   = ensemble_D_mean,
    Pred_EnsStd    = ensemble_D_std,
    Pred_Lower95   = ensemble_D_lower,
    Pred_Upper95   = ensemble_D_upper,
    Pred_Q025      = ensemble_D_q025,
    Pred_Q975      = ensemble_D_q975
)
CSV.write(SAVE_PATH * "weekly_cases_ensemble.csv", df_weekly)
println("✅ weekly_cases_ensemble.csv  ($(nrow(df_weekly)) 行)")

# ============================================================
# 10. 可视化
# ============================================================

# —— Fig 1: 每周新增病例 ——
plt1 = plot(1:n_weeks, Float64.(trainingdata_2);
    label="Observed", lw=2, color=:blue,
    marker=:circle, markersize=2, alpha=0.8,
    xlabel="Week", ylabel="New Cases",
    title="Weekly New Cases — Best Model + Ensemble CI  (kept $n_final/$N_ENSEMBLE)",
    size=(1100, 500), legend=:topright)

plot!(plt1, 1:n_weeks, ensemble_D_q025;
    fillrange=ensemble_D_q975, fillalpha=0.22, fillcolor=:orange,
    label="95% Ensemble CI", linealpha=0)

for col in 1:n_final
    plot!(plt1, 1:n_weeks, all_pred_D[:, col];
        label=(col==1 ? "Ensemble Members" : ""),
        lw=0.7, color=:gray, alpha=0.4)
end

plot!(plt1, 1:n_weeks, Float64.(best_D);
    label="Best Model (#$best_idx)", lw=2.5, color=:red, linestyle=:dash)
plot!(plt1, 1:n_weeks, ensemble_D_mean;
    label="Ensemble Mean", lw=2, color=:darkorange, linestyle=:dot)

display(plt1)
savefig(plt1, SAVE_PATH * "weekly_cases_ensemble_CI.png")

# —— Fig 2: 每日 Beta ——
plt2 = plot(1:n_days, best_beta;
    label="Best Model β (#$best_idx)", lw=2, color=:green,
    xlabel="Day", ylabel="β",
    title="β(t) — Best Model + Ensemble CI  (kept $n_final/$N_ENSEMBLE)",
    size=(1200, 450), legend=:topright)

plot!(plt2, 1:n_days, ensemble_beta_q025;
    fillrange=ensemble_beta_q975, fillalpha=0.22, fillcolor=:green,
    label="95% Ensemble CI", linealpha=0)

for col in 1:n_final
    plot!(plt2, 1:n_days, all_betas[:, col];
        label=(col==1 ? "Ensemble Members" : ""),
        lw=0.5, color=:gray, alpha=0.35)
end

plot!(plt2, 1:n_days, ensemble_beta_mean;
    label="Ensemble Mean", lw=1.5, color=:darkgreen, linestyle=:dot)

for yr in 1:6
    vline!(plt2, [yr*365]; label="", color=:gray, linestyle=:dash, alpha=0.5)
end

display(plt2)
savefig(plt2, SAVE_PATH * "daily_beta_ensemble_CI.png")

# —— Fig 3: 残差 ——
plt3 = plot(layout=(1,2), size=(1300, 450))

residuals_best = Float64.(trainingdata_2) .- Float64.(best_D)
plot!(plt3[1], 1:n_weeks, residuals_best;
    label="Best Model Residuals", lw=1.2, color=:darkblue,
    xlabel="Week", ylabel="Residuals", title="Residuals (Best Model)")
hline!(plt3[1], [0]; label="", color=:gray, linestyle=:dash)

histogram!(plt3[2], residuals_best;
    bins=30, label="Residual Distribution", color=:steelblue, alpha=0.7,
    xlabel="Residuals", ylabel="Frequency", title="Residual Histogram")

display(plt3)
savefig(plt3, SAVE_PATH * "residuals_ensemble.png")

# —— Fig 4: 累积病例 ——
plt4 = plot(1:n_weeks, Float64.(trainingdata_1);
    label="Observed Cumulative", lw=2, color=:blue,
    xlabel="Week", ylabel="Cumulative Cases",
    title="Cumulative Cases (Best Model)", size=(1000, 400))
plot!(plt4, 1:n_weeks, Float64.(best_C);
    label="Best Model", lw=2, color=:red, linestyle=:dash)
display(plt4)
savefig(plt4, SAVE_PATH * "cumulative_ensemble.png")

# ============================================================
# 11. 拟合指标
# ============================================================
println("\n" * "="^60)
println("📊 Goodness of Fit Metrics")
println("-"^40)

obs = Float64.(trainingdata_2)

# 基准模型指标
println("\n【最佳基准模型 #$best_idx】")
println("  NLL  = $(round(best_member.nll, digits=2))")
println("  r    = $(round(best_member.r, digits=3))")
mse_b  = mean((obs .- Float64.(best_D)).^2)
rmse_b = sqrt(mse_b)
mae_b  = mean(abs.(obs .- Float64.(best_D)))
mape_b = mean(abs.(obs .- Float64.(best_D)) ./ max.(obs, 1.0)) * 100
println("  MSE  = $(round(mse_b,  digits=2))")
println("  RMSE = $(round(rmse_b, digits=2))")
println("  MAE  = $(round(mae_b,  digits=2))")
println("  MAPE = $(round(mape_b, digits=2))%")

# 集成不确定性覆盖率
println("\n【Ensemble 不确定性 ($n_final 个成员)】")
in_q  = sum((obs .>= ensemble_D_q025) .& (obs .<= ensemble_D_q975))
in_s  = sum((obs .>= ensemble_D_lower) .& (obs .<= ensemble_D_upper))
println("  95% CI Coverage (quantile): $(round(in_q/n_weeks*100, digits=1))%")
println("  95% CI Coverage (±1.96σ):   $(round(in_s/n_weeks*100, digits=1))%")
avg_width = mean(ensemble_D_q975 .- ensemble_D_q025)
println("  Average CI Width (quantile): $(round(avg_width, digits=2))")

# 集成均值指标
println("\n【Ensemble 均值】")
mse_e  = mean((obs .- ensemble_D_mean).^2)
rmse_e = sqrt(mse_e)
mae_e  = mean(abs.(obs .- ensemble_D_mean))
mape_e = mean(abs.(obs .- ensemble_D_mean) ./ max.(obs, 1.0)) * 100
println("  MSE  = $(round(mse_e,  digits=2))")
println("  RMSE = $(round(rmse_e, digits=2))")
println("  MAE  = $(round(mae_e,  digits=2))")
println("  MAPE = $(round(mape_e, digits=2))%")
println("="^60)