function calculateCS(
        datalist::datalist_struct,
        price::Vector{Float64},
        theta1::Vector{Float64},
        theta2::Vector{Float64},
        Xi::Vector{Float64},
        HH::Int64
    )
    
    X1_new = datalist.X1[:, :];
    X2_new = reshape(price, (:, 1));
    X1_new[:, 2] .= price;
    
    delta = (X1_new * theta1) .+ Xi;
    
    # elasticity
    mu = X2_new * Diagonal(theta2) * datalist.randomDrawMat;
    
    V = delta .+ mu;
    exp_V = exp.(V);
    
    numerator = log.(vec(sum(exp_V, dims = 1)) .+ 1.0);
    
    draw_for_price = datalist.randomDrawMat[1,:];
    alpha_i = - (theta1[2] .+ theta2[1] .* draw_for_price);
    
    CS = mean(numerator ./ alpha_i) .* HH;

    return CS
    
end
