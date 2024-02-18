function generateMileageTransition(
    kappa::Vector{Float64},
    num_mileage_states::Int
    )
    kappa_1 = kappa[1];
    kappa_2 = kappa[2];

    mileage_trans_mat = zeros(num_mileage_states, num_mileage_states, 2);
    for i in 1:num_mileage_states, j in 1:num_mileage_states
        if (i == j)
            mileage_trans_mat[i, j, 1] = 1 - kappa_1 - kappa_2;
        elseif (i == j - 1)
            mileage_trans_mat[i, j, 1] = kappa_1;
        elseif (i == j - 2)
            mileage_trans_mat[i, j, 1] = kappa_2;
        end
    end
    mileage_trans_mat[num_mileage_states - 1, num_mileage_states, 1] = kappa_1 + kappa_2;
    mileage_trans_mat[num_mileage_states, num_mileage_states, 1] = 1;

    mileage_trans_mat[:, :, 2] = repeat(mileage_trans_mat[1, :, 1]', num_mileage_states);

    return mileage_trans_mat
end

