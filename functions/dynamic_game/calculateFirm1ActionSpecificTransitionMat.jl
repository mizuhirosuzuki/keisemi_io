function calculateFirm1ActionSpecificTransitionMat!(
    economyTransitionMat, 
    CCP1Mat,
    stateTransitionMatFirm1ActionSpecific,
    possibleActionArray
    )

    TempMat0 = kron(economyTransitionMat, ones(4, 4));

    TempMat1 = kron(
        ones(1, 2),
        vcat([
            kron(
                CCP1Mat[i, possibleActionArray[i, :, 1] .== 1]',
                [1 1]
            )
            for i in 1:8
        ]...)
    )

    MatAdjustMinus1 = repeat(vcat(zeros(1, 8), ones(1, 8)), 4);
    MatAdjustMinus2 = repeat(vcat(ones(1, 8), zeros(1, 8)), 4)';
    stateTransitionMatFirm1ActionSpecific[:, :, 1] .= (
        TempMat0 .* TempMat1 .* MatAdjustMinus1 .* MatAdjustMinus2
    );


    MatAdjustZero = repeat(vcat(
        repeat([1, 0], 4)',
        repeat([0, 1], 4)',
    ), 4);
    stateTransitionMatFirm1ActionSpecific[:, :, 2] .= (
        TempMat0 .* TempMat1 .* MatAdjustZero
    );


    MatAdjustPlus1 = repeat(vcat(ones(1, 8), zeros(1, 8)), 4);
    MatAdjustPlus2 = repeat(vcat(zeros(1, 8), ones(1, 8)), 4)';
    stateTransitionMatFirm1ActionSpecific[:, :, 3] .= (
        TempMat0 .* TempMat1 .* MatAdjustPlus1 .* MatAdjustPlus2
    );

end
