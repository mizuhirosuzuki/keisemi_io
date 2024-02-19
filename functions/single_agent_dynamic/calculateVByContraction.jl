function calculateVByContraction(
    theta::Vector{Float64},
    beta::Float64,
    trans_mat::Array{Float64, 3},
    states_matrix::Matrix{Int};
    num_states::Int = num_states,
    num_choice::Int = num_choice,
    Euler_const::Float64 = Euler_const
    )

    U = calculateFlowUtil(theta, states_matrix);
    V_old = zeros(num_states);
    V_new = similar(V_old);

    diff = 1000;
    tol_level = 1e-10;

    while (diff > tol_level)
        V_new .= log.(sum(
          exp.(U + beta .* hcat(
            view(trans_mat, :, :, 1) * V_old, 
            view(trans_mat, :, :, 2) * V_old
            )),
          dims = 2)) .+ Euler_const;

        diff = sum(abs.(V_new - V_old));

        V_old .= V_new[:];
    end

    return V_old

end

function calculateVByContraction(
    theta::AbstractArray,
    beta::Float64,
    trans_mat::Array{Float64, 3},
    states_matrix::Matrix{Int};
    num_states::Int = num_states,
    num_choice::Int = num_choice,
    Euler_const::Float64 = Euler_const
    )

    U = calculateFlowUtil(theta, states_matrix);
    V_old = zeros(num_states);

    diff = 1000;
    tol_level = 1e-10;

    while (diff > tol_level)
        V_new = log.(sum(
          exp.(U + beta .* hcat(
            view(trans_mat, :, :, 1) * V_old, 
            view(trans_mat, :, :, 2) * V_old
            )),
          dims = 2)) .+ Euler_const;

        diff = sum(abs.(V_new - V_old));

        V_old = V_new[:];
    end

    return V_old

end

