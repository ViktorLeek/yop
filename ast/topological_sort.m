function [topsort, visited, n_elem] = topological_sort(obj, topsort, visited, n_elem)
% Topological sort of expression graph by a dfs.

% Initialize if second and third args are empty
if nargin == 1
    % topsort = {};
    visited = [];
    topsort = cell(1e4, 1);
    n_elem = 0;
end

end