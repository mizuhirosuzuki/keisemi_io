function updateCCP(
    parameterMat,
    CCP1Mat,
    CCP2Mat,
    TransitionMat,
    global_param
)

    profitFirm1 = calculateFirmProfit(
        parameterMat,
        global_param.entryMat,
        global_param.numFirmsVec,
        global_param.economyVec,
        global_param.possibleActionArray,
        1
    );
    profitFirm2 = calculateFirmProfit(
        parameterMat,
        global_param.entryMat,
        global_param.numFirmsVec,
        global_param.economyVec,
        global_param.possibleActionArray,
        2
    );

    stateTransitionMat = generateStateTransitionMat(
        TransitionMat, 
        CCP1Mat, 
        CCP2Mat, 
        global_param.possibleActionArray
        );

    equilibriumProfitFirm1 = calculateProfitGivenOpponentCCP(
        profitFirm1, 
        CCP2Mat
        );
    equilibriumProfitFirm2 = calculateProfitGivenOpponentCCP(
        profitFirm2, 
        CCP1Mat
        );

    expectedShockUnderBestActionFirm1 = (
        global_param.eulergamma .- replace(log.(CCP1Mat), -Inf => 0)
    );
    expectedShockUnderBestActionFirm2 = (
        global_param.eulergamma .- replace(log.(CCP2Mat), -Inf => 0)
    );

    exanteV1 = vec(
        (I(8) - global_param.beta .* stateTransitionMat) \ 
        sum(
            CCP1Mat .* (
                equilibriumProfitFirm1 + 
                expectedShockUnderBestActionFirm1
                ), dims = 2)
    );

    exanteV2 = vec(
        (I(8) - global_param.beta .* stateTransitionMat) \ 
        sum(
            CCP2Mat .* (
                equilibriumProfitFirm2 + 
                expectedShockUnderBestActionFirm2
                ), dims = 2)
    );

    stateTransitionMatFirm2ActionSpecific = zeros(8, 8, 3);
    stateTransitionMatFirm1ActionSpecific = zeros(8, 8, 3);
    calculateFirm2ActionSpecificTransitionMat!(
        TransitionMat, 
        CCP2Mat,
        stateTransitionMatFirm2ActionSpecific,
        global_param.possibleActionArray
        );
    calculateFirm1ActionSpecificTransitionMat!(
        TransitionMat, 
        CCP1Mat,
        stateTransitionMatFirm1ActionSpecific,
        global_param.possibleActionArray
        );

    updatedCCP1Numerator = zeros(8, 3);
    updatedCCP2Numerator = zeros(8, 3);
    for i in 1:3
        updatedCCP1Numerator[:, i] .= exp.(
            view(equilibriumProfitFirm1, :, i) .+ 
            global_param.beta .* view(stateTransitionMatFirm2ActionSpecific, :, :, i) * 
            exanteV1 .*
            view(global_param.possibleActionArray, :, i, 1)
        );
    end
    CCP1MatUpdated = (
        (
            updatedCCP1Numerator ./ 
            (sum(updatedCCP1Numerator, dims = 2) .- 1)) .* 
        view(global_param.possibleActionArray, :, :, 1)
    );

    for i in 1:3
        updatedCCP2Numerator[:, i] .= exp.(
            view(equilibriumProfitFirm2, :, i) .+ 
            global_param.beta .* view(stateTransitionMatFirm1ActionSpecific, :, :, i) * 
            exanteV2 .*
            view(global_param.possibleActionArray, :, i, 2)
        );
    end
    CCP2MatUpdated = (
        (
            updatedCCP2Numerator ./ 
            (sum(updatedCCP2Numerator, dims = 2) .- 1)
            ) .* 
        view(global_param.possibleActionArray, :, :, 2)
    );

    return [CCP1MatUpdated, CCP2MatUpdated]

end
