function calculateProfitExcludeCompetition(
    param,
    df::DataFrame,
    u_m0::Matrix{Float64},
    u_mIm::Matrix{Float64}
)

    alpha = param[1:8];
    beta = param[9:12];
    delta = param[13];
    rho = param[14];
    
    variableProfit = Matrix(df[:, [:Const, :Menseki, :LogPop, :LogIncome]]) * beta .+ rho .* u_m0;
    fixedCost = Matrix(df[:, [
        :Kyukyu, :Sien, :Hyoka, :DepNeurology, :DepNeurosurgery, :LogNumBeds, :ZeroBedDum, :DaigakuDum
        ]]) * alpha .- (sqrt(1.0 - rho^2) .* u_mIm);
    profitExcludeCompetition = variableProfit - fixedCost;

    return profitExcludeCompetition
end
