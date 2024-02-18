function calculateFlowUtil(
    theta,
    states_matrix::Matrix{Int}
    )

    theta_c = theta[1];
    theta_p = theta[2];

    return hcat(
        - theta_c .* states_matrix[:, 2],
        - theta_p .* states_matrix[:, 1]
    );
end
