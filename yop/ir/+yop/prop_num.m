function value = prop_num(expr)
% PROP_NUM - Propagate numerical values

vars = yop.get_vars(expr);
for k=1:length(vars)   
    vars{k}.m_value = nan(size(vars{k}));
end
value = propagate_value(expr);

end