function addIVColumnsForNestedLogit!(
    data::DataFrame, 
    variable::String
    )
    data[!, Symbol("iv_BLP_own_", variable, "_nest")] = (
        data[:, Symbol(variable, "_sum_own")] - data[:, variable]
    );
    data[!, Symbol("iv_BLP_other_", variable, "_nest")] = (
        data[:, Symbol(variable, "_sum_mkt")] - data[:, Symbol(variable, "_sum_own")]
    );
    data[!, Symbol("iv_GH_own_", variable, "_nest")] = (
        (data[:, :group_n] .- 1) .* data[:, variable].^2 .+ 
        (data[:, Symbol(variable, "_sqr_sum_own")] .- data[:, variable].^2) .- 
        2 .* data[:, variable] .* 
        (data[:, Symbol(variable, "_sum_own")] .- data[:, variable])
    );
    data[!, Symbol("iv_GH_other_", variable, "_nest")] = (
        (data[:, :mkt_n] .- data[:, :group_n]) .* data[:, variable].^2 .+ 
        (data[:, Symbol(variable, "_sqr_sum_mkt")] .- 
        data[:, Symbol(variable, "_sqr_sum_own")]) .- 
        2 .* data[:, variable] .* 
        (data[:, Symbol(variable, "_sum_mkt")] .- 
        data[:, Symbol(variable, "_sum_own")])
    );
end

