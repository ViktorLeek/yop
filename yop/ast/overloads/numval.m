function val = numval(obj)
if isnumeric(obj)
    val = obj;
else
    val = nan(size(obj));
end
end