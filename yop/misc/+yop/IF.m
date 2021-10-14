function val = IF(bool, tval, fval)
% Convenience function. If tval and fval are expensive to compute this
% should be avoided, as both branches might be evaluated. To avoid
% evaluation the argument can be passed as an anonymous function, this can
% be expensive in itself, but can avoid evaluation that leads to errors. An
% empty anonymous function is declared as 'name = @() expr;' it is 
% evaluated by name(). If name is not a function, no harm is done, as '()'
% should not carry side effects in that case.
if bool
    val = tval();
else
    val  = fval();
end
end