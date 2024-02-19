function Estimation_forward_PSD(
    param,
    W1star,
    W2star,
    TransitionMat,
    CCP1Mat,
    CCP2Mat,
    global_param
)

    parameterMat = hcat(
        [param[1]; param[3:4]; [0]; param[5]],
        [param[2]; param[3:4]; [0]; param[5]], 
    );

    NumMaxchoices = 3;

    exanteV1 = W1star' * [parameterMat[:, 1]; 1];
    exanteV2 = W2star' * [parameterMat[:, 2]; 1];

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

    equilibriumProfitFirm1 = calculateProfitGivenOpponentCCP(
        profitFirm1, 
        CCP2Mat
        );
    equilibriumProfitFirm2 = calculateProfitGivenOpponentCCP(
        profitFirm2, 
        CCP1Mat
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

    return sqL2dist(CCP1MatUpdated[:, 2], CCP1Mat[:, 2]) + sqL2dist(CCP2MatUpdated[:, 2], CCP2Mat[:, 2])

end
