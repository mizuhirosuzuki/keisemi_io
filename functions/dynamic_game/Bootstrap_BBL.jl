function Bootstrap_BBL(
    FakeData,
    global_param,
    EVrandom,
    UNIrandom,
    InitialState,
    NumSimMarkets,
    NumSimulations,
    NumSimPeriods,
    NumPerturbations
)

    EstimatedCCP1Mat, EstimatedCCP2Mat = estimateCCPMat(
        FakeData, global_param
    );

    EstimatedTransition = estimateTransition(FakeData);

    output = VSigmaGeneration(
        EstimatedCCP1Mat[:, 2],
        EstimatedCCP2Mat[:, 2],
        EstimatedTransition,
        EVrandom,
        UNIrandom,
        InitialState,
        NumSimMarkets,
        NumSimulations,
        NumSimPeriods,
        global_param
    );

    W1star = output[1];
    W2star = output[2];

    PerturbedCCP1 = reshape(
        clamp.(
            repeat(EstimatedCCP1Mat[:, 2], NumPerturbations) + 
            rand(Normal(0, .1), 8 * NumPerturbations),
            0.001, 0.999
        ),
        (8, NumPerturbations)
    );

    PerturbedCCP2 = reshape(
        clamp.(
            repeat(EstimatedCCP2Mat[:, 2], NumPerturbations) + 
            rand(Normal(0, .1), 8 * NumPerturbations),
            0.001, 0.999
        ),
        (8, NumPerturbations)
    );

    W1_all = zeros(6, NumSimMarkets, NumPerturbations);
    W2_all = zeros(6, NumSimMarkets, NumPerturbations);

    for per in 1:NumPerturbations

        W1_p = VSigmaGeneration(
            PerturbedCCP1[:, per],
            EstimatedCCP2Mat[:, 2],
            EstimatedTransition,
            EVrandom,
            UNIrandom,
            InitialState,
            NumSimMarkets,
            NumSimulations,
            NumSimPeriods,
            global_param
        );
        W1_all[:, :, per] = W1_p[1];

        W2_p = VSigmaGeneration(
            EstimatedCCP1Mat[:, 2],
            PerturbedCCP2[:, per],
            EstimatedTransition,
            EVrandom,
            UNIrandom,
            InitialState,
            NumSimMarkets,
            NumSimulations,
            NumSimPeriods,
            global_param
        );
        W2_all[:, :, per] = W2_p[2];

    end

    initial = [0.3, 0.2, -0.27, 0.45, -2.1];

    result = optimize(
        x -> BBLobjective(
            x,
            NumPerturbations,
            W1star,
            W2star,
            W1_all,
            W2_all
        ),
        initial,
        Optim.Options(show_trace = false)
    )

    return [
        result,
        EstimatedCCP1Mat,
        EstimatedCCP2Mat,
        EstimatedTransition
    ]
end
