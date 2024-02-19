function calculateLikelihoodFromCCP(
  theta,
  CCP::Matrix{Float64},
  data_gen::DataFrame,
  beta::Float64,
  trans_mat::Array{Float64, 3},
  states_matrix::Matrix{Int},
  calculateCCP
  )

  CCP_est = calculateCCP(theta, CCP, beta, trans_mat, states_matrix);

  return sum(log.([CCP_est[x.state, x.action + 1] for x in eachrow(data_gen)]))

end
