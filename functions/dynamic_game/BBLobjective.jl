function BBLobjective(
    theta,
    NumPerturbations,
    W1star,
    W2star,
    W1_all,
    W2_all
    )

    parameterMat = hcat(
        [theta[1]; theta[3:4]; [0]; theta[5]],
        [theta[2]; theta[3:4]; [0]; theta[5]], 
    );

    objvalue = 0;

    for per in 1:NumPerturbations

        temp1 = min.(
            W1star' * [parameterMat[:, 1]; 1] - W1_all[:, :, per]' * [parameterMat[:, 1]; 1], 
            0
            );
        temp2 = min.(
            W2star' * [parameterMat[:, 2]; 1] - W2_all[:, :, per]' * [parameterMat[:, 2]; 1], 
            0
            );

        objvalue += temp1' * temp1 + temp2' * temp2;
    end

    return objvalue
end
