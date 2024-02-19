function calculateRevenue(
        price_cand,
        data,
        datalist,
        delta,
        theta1,
        theta2,
        option
    )

    mc_betado = 3.198 * (1.0 - 1.0 / abs(-2.16720791));

    tempprice = data.price[:];
    tempprice[(data.NameID .== 197) .& (data.year .== 2016)] .= price_cand;
    
    X1_new = datalist.X1[:, :];
    X2_new = datalist.X2[:, :];
    X1_new[:, 2] = tempprice;
    X2_new[:, 1] = tempprice;
    
    org_xi = delta .- datalist.X1 * theta1;
    new_delta = X1_new * theta1 .+ org_xi;
    
    datalist_temp = datalist_struct(
        X1_new, X2_new, datalist.Z, data.share, datalist.marketIndex, 
        data.logit_share, datalist.randomDrawMat, datalist.marketIndexMat
        );
    
    mktshare = calculateMarketShare(theta2, datalist_temp, new_delta);
    
    quant = mktshare .* data.HH;
    revenue = tempprice .* quant;
        
    revenuevec  = revenue[(data.NameID .== 197) .& (data.year .== 2016)];
    revenuevec2 = sum(revenue[in(NIPPYOautoIDvec).(data[:, :NameID]) .& (data.year .== 2016)]);
    
    pivec  = revenuevec  .- mc_betado .* quant[(data.NameID .== 197) .& (data.year .== 2016)];
    pivec2 = revenuevec2 .- mc_betado .* quant[(data.NameID .== 197) .& (data.year .== 2016)];

    if option == "own"
        return(revenuevec[1])
    elseif option == "total"
        return(revenuevec2[1])
    elseif option == "ownpi"
        return(pivec[1])
    elseif option == "totalpi"
        return(pivec2[1])
    end
    
end
