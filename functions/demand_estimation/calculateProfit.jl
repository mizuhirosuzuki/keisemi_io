function calculateProfit(
        Maker::AbstractVector, 
        price::Vector{Float64}, 
        mc::Vector{Float64}, 
        share::Vector{Float64}, 
        HH::Vector{Int64}
    )
    
    dt = DataFrame(
        Maker = Maker,
        profit = (price - mc) .* share .* HH,
        revenue = price.* share .* HH
    )
    
    return combine(
        groupby(dt, :Maker), 
        [:profit, :revenue] .=> sum .=> [:profit, :revenue]
        )
end
