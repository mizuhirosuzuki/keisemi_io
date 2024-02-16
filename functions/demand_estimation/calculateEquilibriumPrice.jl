function calculateEquilibriumPrice(
        datalist::datalist_struct,
        p_ini::Vector{Float64},
        Ownership::BitMatrix,
        theta1::Vector{Float64},
        theta2::Vector{Float64},
        mc::Vector{Float64},
        Xi::Vector{Float64}
    )
    
    lambda = 1e-6;
    p_old = p_ini;
    distance = 10000;

    while (distance > lambda)
        p_new = updatePrice(datalist, p_old, Ownership, theta1, theta2, mc, Xi);
        distance = maximum(abs.(p_new - p_old));
        p_old = p_new[:];
    end

    return p_old
end
