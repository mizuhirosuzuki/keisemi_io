function calculateShareOwnMRI(df)

    output = combine(df, nrow => :Total, :MRIOwnDum => sum => :MRIHos)
    output[!, :PerMRI] = round.(output.MRIHos ./ output.Total .* 100.0, digits = 2)
    
    return Matrix(output)
    
end
