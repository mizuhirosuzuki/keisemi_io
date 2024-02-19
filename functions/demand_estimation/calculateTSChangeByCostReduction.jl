function calculateTSChangeByCostReduction(
        costReductions::Float64,
        costReducingFirms::AbstractVector,
        Ownership::BitMatrix,
        data::DataFrame,
        mc::Vector{Float64},
        datalist::datalist_struct,
        theta1::Vector{Float64},
        theta2::Vector{Float64},
        HH::Vector{Int64},
        p_pre::Vector{Float64},
        pro_rev_pre::DataFrame,
        CS_pre::Float64,
        Xi::Vector{Float64}
    )
   
    mc_new = mc[:];
    mc_new[in(costReducingFirms).(data.Maker)] .*= costReductions;
    
    p_post = calculateEquilibriumPrice(
        datalist, 
        p_pre, 
        Ownership, 
        theta1, 
        theta2, 
        mc_new, 
        Xi
        );
    
    CV = calculateCS(datalist, p_post, theta1, theta2, Xi, HH[1]) - CS_pre;
    
    share_post = simulateMarketShare(datalist, p_post, theta1, theta2, Xi);
    pro_rev_post = calculateProfit(data.Maker, p_post, mc, share_post, HH);
    
    TS_change = CV + sum(pro_rev_post.profit .- pro_rev_pre.profit);
    return TS_change
    
end
