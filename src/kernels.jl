export M4, M5, M6

const σ4 = (2/3, 10/7/π, 1/π)
const σ5 = (1/24, 96/1199/π, 1/20/π)
const σ6 = (1/120, 7/478/π, 1/120/π)


const test_p = p"2x"

const _M6 = @piecewise_polynomial begin
    0.0 => p"0.0"
    1.0 => p"(3 - q)^5 - 6*(2 - q)^5 + 15*(1 - q)^5"
    2.0 => p"(3 - q)^5 - 6*(2 - q)^5"
    3.0 => p"(3 - q)^5"
    _ => p"0.0"
end

#const M6 = _M6 .* σ6