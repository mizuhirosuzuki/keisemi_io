function calculateBerryObjectiveAtCityLevel(
        param,
        df::DataFrame,
        u_m0::Matrix{Float64},
        u_mIm::Matrix{Float64},
        numEntryObs::Vector{Int64}
    )

    profitExcludeCompetition = calculateProfitExcludeCompetition(
        param,
        df,
        u_m0,
        u_mIm
    )

    each_entry_mat = calculateEquilibriumNumEntry(
        df,
        param[13],
        profitExcludeCompetition
    );

    return mean((numEntryObs .- mean(each_entry_mat, dims = 2)).^2)

end

