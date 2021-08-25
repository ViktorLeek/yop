function fn_hdl = get_constructor(obj)
eval(['fn_hdl = @(lhs, rhs) ', class(obj), '(lhs, rhs);']);
end