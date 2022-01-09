function bool = isa_numeric(obj)
if isa(obj, 'function_handle')
    bool = false;
else
    bool = ~isnan(obj);
end
end