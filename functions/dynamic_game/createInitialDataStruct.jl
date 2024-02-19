function createInitialDataStruct(
    profitFirm1,
    profitFirm2,
    TransitionMat,
    global_param
)

    CCP1Stay = repeat([0.5], 8);
    CCP2Stay = repeat([0.5], 8);

    CCP1Mat = (
        hcat(
            1 .- CCP1Stay,
            CCP1Stay,
            1 .- CCP1Stay,
        ) .* 
        global_param.possibleActionArray[:, :, 1]
        );
    CCP2Mat = (
        hcat(
            1 .- CCP2Stay,
            CCP2Stay,
            1 .- CCP2Stay,
        ) .* 
        global_param.possibleActionArray[:, :, 2]
        );

    stateTransitionMat = generateStateTransitionMat(
        TransitionMat, 
        CCP1Mat, 
        CCP2Mat, 
        global_param.possibleActionArray
        );

    equilibriumProfitFirm1 = calculateProfitGivenOpponentCCP(profitFirm1, CCP2Mat);
    equilibriumProfitFirm2 = calculateProfitGivenOpponentCCP(profitFirm2, CCP1Mat);

    expectedShockUnderBestActionFirm1 = (
        global_param.eulergamma .- replace(log.(CCP1Mat), -Inf => 0)
    )
    expectedShockUnderBestActionFirm2 = (
        global_param.eulergamma .- replace(log.(CCP2Mat), -Inf => 0)
    )

    exanteV1 = vec(
        (I(8) - global_param.beta .* stateTransitionMat) \ 
        sum(CCP1Mat .* (equilibriumProfitFirm1 + expectedShockUnderBestActionFirm1), dims = 2)
    );

    exanteV2 = vec(
        (I(8) - global_param.beta .* stateTransitionMat) \ 
        sum(CCP2Mat .* (equilibriumProfitFirm2 + expectedShockUnderBestActionFirm2), dims = 2)
    );


    return array_struct(
        CCP1Mat,
        CCP2Mat,
        stateTransitionMat,
        equilibriumProfitFirm1,
        equilibriumProfitFirm2,
        expectedShockUnderBestActionFirm1,
        expectedShockUnderBestActionFirm2,
        exanteV1,
        exanteV2,
        zeros(8, 3),
        zeros(8, 3),
        zeros(8, 3),
        zeros(8, 3),
        zeros(8, 8, 3),
        zeros(8, 8, 3),
        zeros(8),
        zeros(8)
    );

end
