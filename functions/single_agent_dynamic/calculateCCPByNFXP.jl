function calculateCCPByNFXP(
    theta,
    beta::Float64,
    trans_mat::Array{Float64, 3},
    states_matrix::Matrix{Int};
    num_states::Int = num_states
)

    U = calculateFlowUtil(theta, states_matrix);

    V = calculateVByContraction(theta, beta, trans_mat, states_matrix; num_states = num_states);
    
    CV = U + beta .* hcat(
      view(trans_mat, :, :, 1) * V, 
      view(trans_mat, :, :, 2) * V
      );

    return exp.(CV) ./ sum(exp.(CV), dims = 2);

end
