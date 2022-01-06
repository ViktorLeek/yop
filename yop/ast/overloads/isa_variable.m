function [bool, id, type] = isa_variable(obj)
bool = false(size(obj));
id = zeros(size(obj));
type = yop.var_type.not_var*ones(size(obj));
end