function [type, id] = Type(obj)
type = yop.var_type.not_var*ones(size(obj));
id = zeros(size(obj));
end