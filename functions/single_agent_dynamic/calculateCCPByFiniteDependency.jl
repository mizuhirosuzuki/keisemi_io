function calculateCCPByFiniteDependency(
  theta,
  CCP::Matrix{Float64},
  beta::Float64,
  trans_mat::Array{Float64, 3},
  states_matrix::Matrix{Int}
  )

  U = calculateFlowUtil(theta, states_matrix);

  CV_dif = (
    view(U, :, 2) - view(U, :, 1) +
    beta .* (
      view(trans_mat, :, :, 2) * (- log.(view(CCP, :, 2))) -
      view(trans_mat, :, :, 1) * (- log.(view(CCP, :, 2)))
      )
  );

  prob_buy = exp.(CV_dif) ./ (1 .+ exp.(CV_dif));

  return hcat(1 .- prob_buy, prob_buy)

end
