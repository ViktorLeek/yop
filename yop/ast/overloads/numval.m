function val = numval(obj)
sz = size(obj);
if isnumeric(obj)
    val = obj;
else
    val = nan(sz);
end
end