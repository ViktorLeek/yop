function boolv = isa_numeric(obj)
if isnumeric(obj)
    boolv = true(size(obj));
else
    boolv = false(size(obj));
end
end