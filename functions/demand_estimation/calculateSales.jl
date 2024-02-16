function calculateSales(
        price::Float64,
        year::Int, 
        targetNameID::Int,
        data::DataFrame,
        estimationResult
    )
    
    subData = data[data.year .== year, :];
    subData[(subData[:, :NameID] .== targetNameID), :price] .= price;
    subData[!, :delta] = predict(
        estimationResult,
        subData
    ) + subData.xi_fit;
    subData[!, :denom] .= 1 .+ sum(exp.(subData[:, :delta]));
    subData[!, :pred_sales] = (
        exp.(subData[:, :delta]) ./ subData[:, :denom] .* subData[:, :HH]
    );
    subData = subData[subData.NameID .== targetNameID, :];
    
    return subData.pred_sales[1]

end;

