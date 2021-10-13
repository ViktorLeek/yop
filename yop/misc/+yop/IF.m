function val = IF(bool, tval, fval)
% Convenience function. If tval and fval are expensive to compute this
% should be avoided, as both branches are evaluated.
if bool
    val = tval;
else
    val  = fval;
end
end