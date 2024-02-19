function updatePrice(
    datalist::datalist_struct,
    p_old::Vector{Float64},
    Ownership::BitMatrix,
    theta1::Vector{Float64},
    theta2::Vector{Float64},
    mc::Vector{Float64},
    Xi::Vector{Float64}
)
    
    X1_new = datalist.X1[:, :];
    X2_new = reshape(p_old, (:, 1));
    X1_new[:, 2] .= p_old;
    
    delta = (X1_new * theta1) .+ Xi;
    datalist_new = datalist_struct(
        X1_new, X2_new, datalist.Z, datalist.ShareVec, datalist.marketIndex, 
        datalist.logitshare, datalist.randomDrawMat, datalist.marketIndexMat
        );
    Sharevec = calculateMarketShare(
        theta2, datalist_new, delta
    );
    
    # elasticity

    elasmat = calculateElasticity(
        p_old,
        datalist_new.X2,
        theta1,
        theta2,
        datalist_new.randomDrawMat,
        delta
    );

    Derivative = - elasmat .* Sharevec' ./ p_old;
    Delta = Derivative .* Ownership;
    p_new = mc .+ (Delta \ Sharevec)

    return p_new
    
end
