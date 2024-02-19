function calculateLogLikelihoodStatic(
    theta,
    states_matrix::Matrix{Int},
    data_gen::DataFrame
)
    U = calculateFlowUtil(theta, states_matrix);
    prob_C_stat = exp.(U) ./ sum(exp.(U), dims = 2);
    return sum(
        log.([
                prob_C_stat[x.state, x.action + 1] 
                for x in eachrow(data_gen)
            ])
            )
end

