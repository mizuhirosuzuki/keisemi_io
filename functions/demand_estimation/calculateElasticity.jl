function calculateElasticity(
    price,
    X2,
    theta1,
    theta2,
    randomDrawMat,
    delta
)

    marketSize = size(X2)[1];

    mu = X2 * Diagonal(theta2) * randomDrawMat;
    delta_mu = delta .+ mu;
    exp_delta_mu = exp.(delta_mu);
    denom_outside = exp(0.0);
    denom = repeat(sum(exp_delta_mu, dims = 1), marketSize) .+ denom_outside;

    s_jt_i = exp_delta_mu ./ denom;
    draw_for_price = randomDrawMat[1,:];
    alpha_i = theta1[2] .+ theta2[1] .* draw_for_price;

    mean_s = mean(s_jt_i, dims = 2);

    elasmat = zeros((marketSize, marketSize));
    tmp_vec = similar(alpha_i);

    for j in 1:marketSize, k in 1:marketSize
        if (k != j)
            tmp_vec .= alpha_i .* view(s_jt_i, j, :) .* view(s_jt_i, k, :);
            elasmat[k, j] = (
                - price[k] / mean_s[j] * mean(tmp_vec)
            );
        elseif (k == j)
            tmp_vec .= alpha_i .* view(s_jt_i, j, :) .* (1.0 .- view(s_jt_i, j, :));
            elasmat[k, j] = (
                price[j] / mean_s[j] * mean(tmp_vec)
            );
        end
    end

    return elasmat

end
