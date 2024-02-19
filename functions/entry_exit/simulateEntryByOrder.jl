function simulateEntryByOrder(
    delta,
    cityIndex::BitMatrix,
    profitExcludeCompetition::Matrix{Float64},
    each_entry_mat::Matrix{Int};
    numPotenHos_max::Int = numPotenHos_max
)

    numCity = size(cityIndex)[2];
    numSim = size(profitExcludeCompetition)[2];

    entryMat = Array{Int}(undef, size(profitExcludeCompetition));

    for city in 1:numCity

        profitExcludeCompetitionCity = view(profitExcludeCompetition, cityIndex[:, city], :);

        entryMatCity = profitExcludeCompetitionCity .- delta .* log.(each_entry_mat[city, :] .+ 1)' .>= 0;
        profitNMatCity = profitExcludeCompetitionCity .- delta .* log.(each_entry_mat[city, :])' .>= 0;

        for simIndex in 1:numSim

            numEntry = sum(entryMatCity[:, simIndex]);
            hospitalIndex = 1;

            while numEntry < each_entry_mat[city, simIndex]
                if (
                    (entryMatCity[hospitalIndex, simIndex] == 0 && profitNMatCity[hospitalIndex, simIndex] == 1)
                )
                    entryMatCity[hospitalIndex, simIndex] = 1;
                    numEntry += 1;
                end
                hospitalIndex += 1;
            end

        end
        entryMat[cityIndex[:, city], :] .= entryMatCity;
    end

    return entryMat
end
