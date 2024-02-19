function generateData(
    consumer_id::Int,
    V_CS::Matrix{Float64},
    trans_mat::Array{Float64, 3},
    price_dist_steady::Vector{Float64};
    num_period::Int = num_period,
)

    state_vec = zeros(Int, num_period);
    action_vec = zeros(Int, num_period);
    period_vec = 1:num_period;

    eps_type1 = reshape(
        rand(GeneralizedExtremeValue(0, 1, 0), num_period * 2),
        (num_period, 2)
    );

    state_vec[1] = sample(ProbabilityWeights(price_dist_steady));

    for t in 1:(num_period - 1)

        state_id_today = state_vec[t];

        if (
            V_CS[:, 1][state_id_today] + eps_type1[t, 1] > 
            V_CS[:, 2][state_id_today] + eps_type1[t, 2]
            )
            action_vec[t] = 0;
            state_vec[t + 1] = sample(
                ProbabilityWeights(trans_mat[state_id_today, :, 1])
                );
        else
            action_vec[t] = 1;
            state_vec[t + 1] = sample(
                ProbabilityWeights(trans_mat[state_id_today, :, 2])
                );
        end

    end

    return DataFrame(
        state = state_vec, 
        action = action_vec, 
        period = period_vec,
        consumer_id = consumer_id
        )
end
