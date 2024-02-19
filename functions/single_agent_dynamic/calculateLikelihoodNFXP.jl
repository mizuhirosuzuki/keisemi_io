function calculateLikelihoodNFXP(
    theta,
    data_gen::DataFrame,
    beta::Float64,
    trans_mat::Array{Float64, 3},
    states_matrix::Matrix{Int}
)
  CCP = calculateCCPByNFXP(theta, beta, trans_mat, states_matrix);

  return sum(log.([CCP[x.state, x.action + 1] for x in eachrow(data_gen)]))

end
