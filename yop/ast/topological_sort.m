function [topsort, n_elem, visited] = ...
    topological_sort(obj, visited, topsort, n_elem)
% Topological sort of expression graph by a dfs.

% if nargin == 1
%     % Start new sort
%     visited = [];
%     topsort = cell( ...
%         yop.constants().topsort_preallocation_size, 1);
%     n_elem = 0;
% end

switch nargin
    case 1
        % Start new sort
        visited = [];
        topsort = cell(yop.constants().topsort_preallocation_size, 1);
        n_elem = 0;
        
    case 2
        % Semi-warm start
        % In semi-warm start 'visited' is already provided, but no elements
        % are sorted. This is for instance useful for finding all variables
        % in a number of expressions that are suspected to contain common
        % subexpressions.
        topsort = cell(yop.constants().topsort_preallocation_size, 1);
        n_elem = 0;
end

n_elem = n_elem + 1;
topsort{n_elem} = obj;

end