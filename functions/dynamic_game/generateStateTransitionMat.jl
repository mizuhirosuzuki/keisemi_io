function generateStateTransitionMat(
    economyTransitionMat, 
    CCP1Mat, 
    CCP2Mat,
    possibleActionArray
    )

    TempMat0 = kron(economyTransitionMat, ones(4, 4));

    TempMat1 = hcat([
        kron(
            vcat([
            CCP1Mat[i, possibleActionArray[i, :, 1] .== 1]'
            for i in 1:8
        ]...),
        [1 1]
    ) for j in 1:2
    ]...);

    TempMat2 = vcat([
        CCP2Mat[i, possibleActionArray[i, :, 2] .== 1]'
        for i in 1:8
    ]...);
    TempMat2 = kron(ones(1, 4), TempMat2);

    return TempMat0 .* TempMat1 .* TempMat2;

end;

