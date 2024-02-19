function calculateBerryObjectiveAtHospitalLevel(
        param,
        df::DataFrame,
        u_m0::Matrix{Float64},
        u_mIm::Matrix{Float64},
        cityIndex::BitMatrix,
        simulateEntry
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

    entryMat = simulateEntry(
        param[13],
        cityIndex,
        profitExcludeCompetition,
        each_entry_mat
    )

    return mean((df.MRIOwnDum .- mean(entryMat, dims = 2)).^2);

end
