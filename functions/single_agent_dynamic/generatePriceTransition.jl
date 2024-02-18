function generatePriceTransition!(
    lambda::Vector{Float64},
    num_price_states::Int,
    price_trans_mat::Array{Float64, 2}
    )

    price_trans_mat[:, :] .= 0;

    price_trans_mat[1, 2:end] = lambda[1:(num_price_states - 1)];
    price_trans_mat[1, 1] = 1 - sum(price_trans_mat[1, :]);

    for i in 2:(num_price_states - 1)
        price_trans_mat[i, 1:(i - 1)] = lambda[
            ((i - 1) * (num_price_states - 1) + 1):((i - 1) * (num_price_states - 1) + (i - 1))
            ];
        price_trans_mat[i, (i + 1):end] = lambda[
            ((i - 1) * (num_price_states - 1) + i):(i * (num_price_states - 1))
            ];
        price_trans_mat[i, i] = 1 - sum(price_trans_mat[i, :]);
    end

    price_trans_mat[num_price_states, 1:(end - 1)] = lambda[((num_price_states - 1) * (num_price_states - 1) + 1):end];
    price_trans_mat[num_price_states, num_price_states] = 1 - sum(price_trans_mat[num_price_states, :]);

end
