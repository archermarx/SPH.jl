export M4, M5, M6

const σ4 = (2/3, 10/7/π, 1/π)
const σ5 = (1/24, 96/1199/π, 1/20/π)
const σ6 = (1/120, 7/478/π, 1/120/π)

const M4 = σ4 .* @piecewise_polynomial begin
    0.0 => p"0.0"
    1.0 => p"0.25*(2-q)^3 - (1 - q)^3"
    2.0 => p"0.25*(2-q)^3"
    _   => p"0.0"
end

const M5 = σ5 .* @piecewise_polynomial begin
    0.0 => p"0.0"
    0.5 => p"(2.5 - q)^4 - 5*(1.5 - q)^4 + 10*(0.5 - q)^4"
    1.5 => p"(2.5 - q)^4 - 5*(1.5 - q)^4"
    2.5 => p"(2.5 - q)^4"
    _   => p"0.0"
end

const M6 = σ6 .* @piecewise_polynomial begin
    0.0 => p"0.0"
    1.0 => p"(3 - q)^5 - 6*(2 - q)^5 + 15*(1 - q)^5"
    2.0 => p"(3 - q)^5 - 6*(2 - q)^5"
    3.0 => p"(3 - q)^5"
    _   => p"0.0"
end