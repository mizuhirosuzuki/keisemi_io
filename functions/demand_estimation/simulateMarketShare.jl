function simulateMarketShare(
        datalist::datalist_struct,
        p::Vector{Float64},
        theta1::Vector{Float64},
        theta2::Vector{Float64},
        Xi::Vector{Float64}
    )
    
    X1_new = datalist.X1[:, :];
    X2_new = reshape(p[:], (:, 1));
    X1_new[:, 2] .= p;
    
    delta = (X1_new * theta1) .+ Xi;
    datalist_new = datalist_struct(
        X1_new, X2_new, datalist.Z, datalist.ShareVec, datalist.marketIndex, 
        datalist.logitshare, datalist.randomDrawMat, datalist.marketIndexMat
        );
    Sharevec = calculateMarketShare(
        theta2, datalist_new, delta
    );
    
    return(Sharevec)
    
end
