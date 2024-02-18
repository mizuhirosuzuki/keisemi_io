function calculateCCPByMatrixInversion(
  theta,
  CCP::Matrix{Float64},
  beta::Float64,
  trans_mat::Array{Float64, 3},
  states_matrix::Matrix{Int};
  num_states::Int = num_states,
  Euler_const::Float64 = Euler_const
)

  U = calculateFlowUtil(theta, states_matrix);

  psi = Euler_const .- log.(CCP);

  V = (
    diagm(ones(num_states)) - 
    beta .* (
      view(CCP, :, 1) .* view(trans_mat, :, :, 1) + 
      view(CCP, :, 2) .* view(trans_mat, :, :, 2)
      )
  ) \ sum(CCP .* (U + psi), dims = 2);

  CV = U + beta .* hcat(
    view(trans_mat, :, :, 1) * V, 
    view(trans_mat, :, :, 2) * V
    );

  return exp.(CV) ./ sum(exp.(CV), dims = 2)

end

