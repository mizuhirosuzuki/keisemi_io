function calculateFirm2ActionSpecificTransitionMat!(
    economyTransitionMat, 
    CCP2Mat,
    stateTransitionMatFirm2ActionSpecific,
    possibleActionArray
    )

    TempMat0 = kron(economyTransitionMat, ones(4, 4));

    TempMat2 = hcat([
        vcat([
            kron(
                [1 1],
                CCP2Mat[i, possibleActionArray[i, :, 2] .== 1]'
            )
            for i in 1:8
        ]...)
        for j in 1:2
    ]...)

    MatAdjustMinus1 = repeat(vcat(
        zeros(2, 8), ones(2, 8)
    ), 2);
    MatAdjustMinus2 = repeat(vcat(
        ones(2, 8), zeros(2, 8)
    ), 2)';
    stateTransitionMatFirm2ActionSpecific[:, :, 1] .= (
        TempMat0 .* TempMat2 .* MatAdjustMinus1 .* MatAdjustMinus2
    );


    ForZero = diagm(ones(2));
    MatAdjustZero = kron(ForZero, ones(2, 2));
    MatAdjustZero = hcat(MatAdjustZero, MatAdjustZero);
    MatAdjustZero = vcat(MatAdjustZero, MatAdjustZero);
    stateTransitionMatFirm2ActionSpecific[:, :, 2] .= (
        TempMat0 .* TempMat2 .* MatAdjustZero;
    );


    MatAdjustPlus1 = repeat(vcat(
        ones(2, 8), zeros(2, 8)
    ), 2);
    MatAdjustPlus2 = repeat(vcat(
        zeros(2, 8), ones(2, 8)
    ), 2)';
    stateTransitionMatFirm2ActionSpecific[:, :, 3] .= (
        TempMat0 .* TempMat2 .* MatAdjustPlus1 .* MatAdjustPlus2
    );

end

