function calculateFirmProfit(
    parameterMat,
    entryMat::Matrix{Int},
    numFirmsVec::Vector{Int},
    economyVec::Vector{Int},
    possibleActionArray::Array{Int, 3},
    firmIndex::Int
    )

    variableProfit = entryMat[:, firmIndex] .* (
        parameterMat[1, firmIndex] .+ (
            parameterMat[2, firmIndex] .* (numFirmsVec .== 2) 
        ) .+ (
            parameterMat[3, firmIndex] .* (economyVec .== 1) 
        )
    )

    fixedCost = repeat(
        [parameterMat[4, firmIndex] 0 parameterMat[5, firmIndex]],
        8
    );

    return (variableProfit .+ fixedCost) .* possibleActionArray[:, :, firmIndex]
end;
