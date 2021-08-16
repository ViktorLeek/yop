function out = tester(tmp)
switch tmp
    case 'yop.ast_le'
        out = 'less than or equal';
    otherwise
        out = 'unknown input';
end
end