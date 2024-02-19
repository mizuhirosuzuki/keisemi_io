function calculateMarketShare(
        theta2,
        datalist::datalist_struct,
        delta
    )
        
    mu = datalist.X2 * Diagonal(theta2) * datalist.randomDrawMat;
    
    delta_mu = delta .+ mu;
    exp_delta_mu = exp.(delta_mu);
    denom_outside = exp(0);
    
    denom_temp = (datalist.marketIndexMat' * exp_delta_mu) .+ denom_outside;
    denom = datalist.marketIndexMat * denom_temp;
    
    s_jt_i = exp_delta_mu ./ denom;
    s_jt = vec(mean(s_jt_i, dims = 2));
    
    return s_jt
    
end