function calculateGMMObjective(
        theta2,
        datalist::datalist_struct,
        delta_ini::Vector{Float64}
    )
    
    delta = calculateMeanUtil(theta2, datalist, delta_ini);
    W = inv(datalist.Z' * datalist.Z);
    
    beta_hat = (
        (datalist.X1' * datalist.Z * W * datalist.Z' * datalist.X1) \ 
        (datalist.X1' * datalist.Z * W * datalist.Z' * delta)
    );
    
    Xi = delta - datalist.X1 * beta_hat;
    
    return Xi' * datalist.Z * W * datalist.Z' * Xi
        
end    

