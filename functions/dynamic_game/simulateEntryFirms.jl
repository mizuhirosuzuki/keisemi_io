function simulateEntryFirms(
    TransitionMat,
    CCP1Mat,
    CCP2Mat,
    global_param,
    NumSimPeriods
)

    stateTransitionMat = generateStateTransitionMat(
        TransitionMat, 
        CCP1Mat, 
        CCP2Mat, 
        global_param.possibleActionArray
        );

    initial_state = [1; zeros(7)];

    transitionpath = zeros(NumSimPeriods, 8);

    for tt in 1:NumSimPeriods
        transitionpath[tt, :] = (
            tt == 1 ? initial_state :
            stateTransitionMat' * transitionpath[tt - 1, :]
        )
    end

    n1 = sum(
        transitionpath .* 
        repeat([0, 0, 1, 1, 0, 0, 1, 1]', NumSimPeriods), 
        dims = 2
        );
    n2 = sum(
        transitionpath .* 
        repeat([0, 1, 0, 1, 0, 1, 0, 1]', NumSimPeriods), 
        dims = 2
        );

    return [n1, n2]

end

