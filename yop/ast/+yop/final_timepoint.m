function tf = final_timepoint(sz)
% Used for propagating timepoints in isa_timepoint queries
if nargin == 0
    sz = 1;
end
tf = inf(sz);
end