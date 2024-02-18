function calculateLogLikelihoodDynamic(
    theta,
    beta::Float64,
    trans_mat::Array{Float64, 3},
    states_matrix::Matrix{Int},
    data_gen::DataFrame
)

    EV = calculateEV(theta, beta, trans_mat, states_matrix);
    
    U = calculateFlowUtil(theta, states_matrix);
    V_CS = U + beta .* EV;

    prob_C = exp.(V_CS) ./ sum(exp.(V_CS), dims = 2);

    return sum(log.([prob_C[x.state, x.action + 1] for x in eachrow(data_gen)]))
end
