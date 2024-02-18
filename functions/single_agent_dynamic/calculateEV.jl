function calculateEV(
    theta::Float64,
    beta::Float64,
    trans_mat::Array{Float64, 3},
    states_matrix::Matrix{Int};
    num_states::Int = num_states,
    num_choice::Int = num_choice,
    Euler_const::Float64 = Euler_const
    )

    U = calculateFlowUtil(theta, states_matrix);

    EV_old = zeros(num_states, num_choice);
    EV_new = similar(EV_old);

    diff = 1000;
    tol_level = 1e-10;

    while (diff > tol_level)
        EV_new[:, 1] .= (
            Euler_const .+ 
            view(trans_mat, :, :, 1) * 
            log.(sum(exp.(U .+ beta .* EV_old), dims = 2))
        )
        EV_new[:, 2] .= (
            Euler_const .+ 
            view(trans_mat, :, :, 2) * 
            log.(sum(exp.(U .+ beta .* EV_old), dims = 2))
        )

        # EV_new = hcat(
        #     Euler_const .+ 
        #     view(trans_mat, :, :, 1) * 
        #     log.(sum(exp.(U .+ beta .* EV_old), dims = 2)),
        #     Euler_const .+ 
        #     view(trans_mat, :, :, 2) * 
        #     log.(sum(exp.(U .+ beta .* EV_old), dims = 2))
        # );

        diff = sum(abs.(EV_new - EV_old));

        EV_old .= EV_new[:, :];
    end

    return EV_old

end


function calculateEV(
    theta::AbstractVector,
    beta::Float64,
    trans_mat::Array{Float64, 3},
    states_matrix::Matrix{Int};
    num_states::Int = num_states,
    num_choice::Int = num_choice,
    Euler_const::Float64 = Euler_const
    )

    U = calculateFlowUtil(theta, states_matrix);

    EV_old = zeros(num_states, num_choice);

    diff = 1000;
    tol_level = 1e-10;

    while (diff > tol_level)

        EV_new = hcat(
            Euler_const .+ 
            view(trans_mat, :, :, 1) * 
            log.(sum(exp.(U .+ beta .* EV_old), dims = 2)),
            Euler_const .+ 
            view(trans_mat, :, :, 2) * 
            log.(sum(exp.(U .+ beta .* EV_old), dims = 2))
        );

        diff = sum(abs.(EV_new - EV_old));

        EV_old = EV_new[:, :];
    end

    return EV_old

end

