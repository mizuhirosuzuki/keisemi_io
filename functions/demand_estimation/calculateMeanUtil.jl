function calculateMeanUtil(
        theta2,
        datalist::datalist_struct,
        delta_ini
    )
    
    tol = 1e-11;
    norm = 1e+10;

    delta_old = delta_ini[:];
    exp_delta_old = exp.(delta_old);

    iter = 0;
        
    while ((norm > tol) & (iter < 1000))
        
        pred_mkt_share = calculateMarketShare(
            theta2, 
            datalist, 
            delta_old
            );
        
        exp_delta = exp_delta_old .* datalist.ShareVec ./ pred_mkt_share;
        
        norm = maximum(abs.(exp_delta .- exp_delta_old));
        
        exp_delta_old = exp_delta[:];
        delta_old = log.(exp_delta_old);
        iter += 1;
        
    end

    return delta_old;
    
end
