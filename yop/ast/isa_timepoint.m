function [bool, tp] = isa_timepoint(obj)
% Convservative implementation. Only that which isa_variable will also
% propagate timepoint information. The reason is that it is mainly used to
% determine whether a box constraint applies to the entire horizon, or 
bool = false(size(obj));
tp = zeros(size(obj));
end