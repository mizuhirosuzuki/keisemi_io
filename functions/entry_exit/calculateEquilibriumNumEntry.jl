function calculateEquilibriumNumEntry(
    df::DataFrame,
    delta::Float64,
    profitExcludeCompetition::Matrix{Float64};
    numPotenHos_max::Int = numPotenHos_max
)

    function each_entry_func(i)
        entry_decision_df = hcat(
            df[:, [:CityCode]], 
            DataFrame((profitExcludeCompetition .- delta .* log(i) .>= 0), :auto)
            );
        return (Matrix(
            combine(
                groupby(entry_decision_df, :CityCode), 
                Not(:CityCode) .=> sum
                )[:, Not(:CityCode)]
                ) .>= i) .* i
    end
    
    each_entry_mat = max.([
        each_entry_func(i) for i in 1:numPotenHos_max
    ]...);

    return each_entry_mat
end
