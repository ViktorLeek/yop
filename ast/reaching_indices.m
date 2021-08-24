function res = reaching_indices(expr)

% Get the variables that make up the expression
vars = get_variables(expr);

% Enumerate all elements
idx0 = 1; % Index to start counting from
var_enum = cell(size(vars)); % indices for the respective variables

for k=1:length(vars)
    % Enumerate variable k
    vk = vars{k};
    sk = size(vk);
    nk = prod(sk);                         % number of elements
    var_enum{k} = idx0:(idx0+n-1);         % enumerate all elements
    vk.m_value = reshape(var_enum{k}, sz); % Necessary to evaluate properly
    idx0 = idx0 + nk;
end

% The value of the node, given the enumeration as input
ridx = forward_evaluate(expr);

% Elements that are not variables are set to some value the enumeration
% does not take.
ridx(~isa_variable(expr)) = -1;

% Separate the indices on a per variable basis
idx_r_var = cell(size(vars));
for k=1:length(vars)
    enum_k = var_enum{k};
    idx_r_var{k} = ridx(ridx >= enum_k(1) & ridx <= enum_k(end));
end


% Compute absolute adress into reaching indices vector
abs_idx = cell(size(vars));
for k=1:length(vars)
    abs_k = false(size(ridx));
    reaching_k = idx_r_var{k};
    for n=1:length(reaching_k)
        abs_k(ridx == reaching_k(n)) = true;
    end
    abs_idx{k} = abs_k;
end

res.ridx = ridx; % reaching indixes
res.vars = vars; % variables
res.vidx = var_enum; % indices for that variable

end