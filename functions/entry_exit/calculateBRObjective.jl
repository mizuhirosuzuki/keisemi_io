function calculateBRObjective(
        params,
        dataset::DataFrame,
        hospitalNumCap::Int
    )

    alpha = [params[1]; -params[2:hospitalNumCap]];
    gamma = params[hospitalNumCap + 1];
    
    NumMRI = dataset.NumMRI;
    numObs = nrow(dataset);
    
    V = LowerTriangular(ones(hospitalNumCap) .* alpha');
        
    VV = (V * ones(hospitalNumCap, numObs))';
    
    profit = dataset.Pop .* VV .- gamma;
    
    phi = cdf.(Normal(0, 1), profit);
    
    mat = hcat(
        1.0 .- phi[:, 1], 
        phi[:, 1:(hospitalNumCap - 1)] .- phi[:, 2:hospitalNumCap], 
        phi[:, hospitalNumCap]
        );
    
    ml = log.(
        reduce(.+, [(dataset.NumMRI .== i) .* mat[:, (i + 1)] for i = 0:hospitalNumCap])
        );
    
    return - sum(ml) / numObs
    
end
