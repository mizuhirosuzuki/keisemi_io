function VSigmaGeneration(
    CCP1,
    CCP2,
    EstimatedTransition,
    EVrandom,
    UNIrandom,
    InitialState,
    NumSimMarkets,
    NumSimulations,
    NumSimPeriods,
    global_param
)

    W1 = zeros(6, NumSimMarkets, NumSimulations);
    W2 = zeros(6, NumSimMarkets, NumSimulations);

    ThresholdValue1 = log.(CCP1) - log.(1 .- CCP1);
    ThresholdValue2 = log.(CCP2) - log.(1 .- CCP2);

    W1Seed = zeros(6);
    W2Seed = zeros(6);

    for mrkt in 1:NumSimMarkets, sim in 1:NumSimulations

        State = InitialState[mrkt, 2];

        for t in 1:NumSimPeriods
        
            DiffEpsilon1inc = -EVrandom[mrkt, t, 1, sim, State, 2] + EVrandom[mrkt, t, 1, sim, State, 3];
            DiffEpsilon1dec = -EVrandom[mrkt, t, 1, sim, State, 2] + EVrandom[mrkt, t, 1, sim, State, 1];

            DiffEpsilon2inc = -EVrandom[mrkt, t, 2, sim, State, 2] + EVrandom[mrkt, t, 2, sim, State, 3];
            DiffEpsilon2dec = -EVrandom[mrkt, t, 2, sim, State, 2] + EVrandom[mrkt, t, 2, sim, State, 1];

            if State in [1, 2, 5, 6]
                a1 = (ThresholdValue1[State] < DiffEpsilon1inc);
                e1 = (
                    (1 - a1) * EVrandom[mrkt, t, 1, sim, State, 2] +
                    a1 * EVrandom[mrkt, t, 1, sim, State, 3]
                );
            else
                a1 = (ThresholdValue1[State] < DiffEpsilon1dec);
                e1 = (
                    (1 - a1) * EVrandom[mrkt, t, 1, sim, State, 2] +
                    a1 * EVrandom[mrkt, t, 1, sim, State, 1]
                );
            end

            if State in [1, 3, 5, 7]
                a2 = (ThresholdValue2[State] < DiffEpsilon2inc);
                e2 = (
                    (1 - a2) * EVrandom[mrkt, t, 2, sim, State, 2] +
                    a2 * EVrandom[mrkt, t, 2, sim, State, 3]
                );
            else
                a2 = (ThresholdValue2[State] < DiffEpsilon2dec);
                e2 = (
                    (1 - a2) * EVrandom[mrkt, t, 2, sim, State, 2] +
                    a2 * EVrandom[mrkt, t, 2, sim, State, 1]
                );
            end

            if State == 1
                W1Seed .= [0, 0, 0, 0, a1, e1];
                W2Seed .= [0, 0, 0, 0, a2, e2];

                NextState = 1 + 2 * a1 + a2;
                TransitionSeed = (UNIrandom[mrkt, t, sim] > EstimatedTransition[1, 1]);
                NextState = 4 * TransitionSeed + NextState;
            elseif State == 2
                W1Seed .= [0, 0, 0, 0, a1, e1];
                W2Seed .= [1, 0, 1, a2, 0, e2];

                NextState = 2 + 2 * a1 - a2;
                TransitionSeed = (UNIrandom[mrkt, t, sim] > EstimatedTransition[1, 1]);
                NextState = 4 * TransitionSeed + NextState;
            elseif State == 3
                W1Seed .= [1, 0, 1, a1, 0, e1];
                W2Seed .= [0, 0, 0, 0, a2, e2];

                NextState = 3 - 2 * a1 + a2;
                TransitionSeed = (UNIrandom[mrkt, t, sim] > EstimatedTransition[1, 1]);
                NextState = 4 * TransitionSeed + NextState;
            elseif State == 4
                W1Seed .= [1, 1, 1, a1, 0, e1];
                W2Seed .= [1, 1, 1, a2, 0, e2];

                NextState = 4 - 2 * a1 - a2;
                TransitionSeed = (UNIrandom[mrkt, t, sim] > EstimatedTransition[1, 1]);
                NextState = 4 * TransitionSeed + NextState;
            elseif State == 5
                W1Seed .= [0, 0, 0, 0, a1, e1];
                W2Seed .= [0, 0, 0, 0, a2, e2];

                NextState = 1 + 2 * a1 + a2;
                TransitionSeed = (UNIrandom[mrkt, t, sim] < EstimatedTransition[2, 2]);
                NextState = 4 * TransitionSeed + NextState;
            elseif State == 6
                W1Seed .= [0, 0, 0, 0, a1, e1];
                W2Seed .= [1, 0, 0, a2, 0, e2];

                NextState = 2 + 2 * a1 - a2;
                TransitionSeed = (UNIrandom[mrkt, t, sim] < EstimatedTransition[2, 2]);
                NextState = 4 * TransitionSeed + NextState;
            elseif State == 7
                W1Seed .= [1, 0, 0, a1, 0, e1];
                W2Seed .= [0, 0, 0, 0, a2, e2];

                NextState = 3 - 2 * a1 + a2;
                TransitionSeed = (UNIrandom[mrkt, t, sim] < EstimatedTransition[2, 2]);
                NextState = 4 * TransitionSeed + NextState;
            elseif State == 8
                W1Seed .= [1, 1, 0, a1, 0, e1];
                W2Seed .= [1, 1, 0, a2, 0, e2];

                NextState = 4 - 2 * a1 - a2;
                TransitionSeed = (UNIrandom[mrkt, t, sim] < EstimatedTransition[2, 2]);
                NextState = 4 * TransitionSeed + NextState;
            end

            W1[:, mrkt, sim] .+= (global_param.beta^(t - 1)).*W1Seed;
            W2[:, mrkt, sim] .+= (global_param.beta^(t - 1)).*W2Seed;

            State = NextState;
        end

    end

    W1out = mean(W1, dims = 3)[:, :, 1];
    W2out = mean(W2, dims = 3)[:, :, 1];

    return [W1out, W2out]

end

