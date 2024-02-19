function calculateMPE!(
    data_struct,
    TransitionMat,
    profitFirm1,
    profitFirm2,
    global_param
)

    diffExanteV = 1;
    iter = 0;

    while (diffExanteV > 1e-12)

        calculateFirm2ActionSpecificTransitionMat!(
            TransitionMat, 
            data_struct.CCP2Mat,
            data_struct.stateTransitionMatFirm2ActionSpecific,
            global_param.possibleActionArray
            );
        calculateFirm1ActionSpecificTransitionMat!(
            TransitionMat, 
            data_struct.CCP1Mat,
            data_struct.stateTransitionMatFirm1ActionSpecific,
            global_param.possibleActionArray
            );

        for i in 1:3
            data_struct.updatedCCP1Numerator[:, i] .= exp.(
                view(data_struct.equilibriumProfitFirm1, :, i) .+ 
                global_param.beta .* view(data_struct.stateTransitionMatFirm2ActionSpecific, :, :, i) * 
                data_struct.exanteV1 .*
                view(global_param.possibleActionArray, :, i, 1)
            );
        end
        data_struct.CCP1MatUpdated .= (
            (
                data_struct.updatedCCP1Numerator ./ 
                (sum(data_struct.updatedCCP1Numerator, dims = 2) .- 1)) .* 
            view(global_param.possibleActionArray, :, :, 1)
        );

        for i in 1:3
            data_struct.updatedCCP2Numerator[:, i] .= exp.(
                view(data_struct.equilibriumProfitFirm2, :, i) .+ 
                global_param.beta .* view(data_struct.stateTransitionMatFirm1ActionSpecific, :, :, i) * 
                data_struct.exanteV2 .*
                view(global_param.possibleActionArray, :, i, 2)
            );
        end
        data_struct.CCP2MatUpdated .= (
            (
                data_struct.updatedCCP2Numerator ./ 
                (sum(data_struct.updatedCCP2Numerator, dims = 2) .- 1)
                ) .* 
            view(global_param.possibleActionArray, :, :, 2)
        );

        data_struct.stateTransitionMat .= generateStateTransitionMat(
            TransitionMat, 
            data_struct.CCP1MatUpdated, 
            data_struct.CCP2MatUpdated, 
            global_param.possibleActionArray
            );

        data_struct.equilibriumProfitFirm1 .= calculateProfitGivenOpponentCCP(
            profitFirm1, 
            data_struct.CCP2MatUpdated
            );
        data_struct.equilibriumProfitFirm2 .= calculateProfitGivenOpponentCCP(
            profitFirm2, 
            data_struct.CCP1MatUpdated
            );

        data_struct.expectedShockUnderBestActionFirm1 .= (
            global_param.eulergamma .- replace(log.(data_struct.CCP1MatUpdated), -Inf => 0)
        )
        data_struct.expectedShockUnderBestActionFirm2 .= (
            global_param.eulergamma .- replace(log.(data_struct.CCP2MatUpdated), -Inf => 0)
        )

        data_struct.exanteV1Updated .= vec(
            (I(8) - global_param.beta .* data_struct.stateTransitionMat) \ 
            sum(
                data_struct.CCP1MatUpdated .* (
                    data_struct.equilibriumProfitFirm1 + 
                    data_struct.expectedShockUnderBestActionFirm1
                    ), dims = 2)
        );

        data_struct.exanteV2Updated .= vec(
            (I(8) - global_param.beta .* data_struct.stateTransitionMat) \ 
            sum(
                data_struct.CCP2MatUpdated .* (
                    data_struct.equilibriumProfitFirm2 + 
                    data_struct.expectedShockUnderBestActionFirm2
                    ), dims = 2)
        );

        diffExanteV = sum(
            (data_struct.exanteV1Updated - data_struct.exanteV1).^2
            ) + sum(
                (data_struct.exanteV2Updated - data_struct.exanteV2).^2
                );

        data_struct.exanteV1 .= data_struct.exanteV1Updated[:];
        data_struct.exanteV2 .= data_struct.exanteV2Updated[:];

        data_struct.CCP1Mat .= data_struct.CCP1MatUpdated[:, :];
        data_struct.CCP2Mat .= data_struct.CCP2MatUpdated[:, :];

        iter += 1;

    end

end

