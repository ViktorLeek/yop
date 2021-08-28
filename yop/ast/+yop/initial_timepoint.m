function t0 = initial_timepoint(sz)
% Used for propagating timepoints in isa_timepoint queries
if nargin == 0
    sz = 1;
end
t0 = -inf(sz);
end