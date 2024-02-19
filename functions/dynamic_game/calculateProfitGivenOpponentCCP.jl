function calculateProfitGivenOpponentCCP(profitMat, CCPOpponent)

    return hcat([
        sum(profitMat[:, i] .* CCPOpponent, dims = 2)
        for i in 1:3
    ]...)

end

