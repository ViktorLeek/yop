function [T_out] = T_exhaust_mix(T_t,T_em,W_t,W_wg)
    T_out = (T_em .* W_wg + T_t .* W_t) ./ (W_t + W_wg);
end