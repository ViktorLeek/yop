function yop_control(varargin)

deg = 0;
n_var = 0;
for k=1:length(varargin)
    if strcmp(varargin{k}, 'deg')
        filt = replace(varargin{k+1}, ...
            {'[', ']', ', ', ' ,', ','}, ...
            {'' ,  '',  ' ',  ' ', ' '} ...
            );
        tmp = textscan(filt, '%f', 'delimiter', ' ');
        deg = tmp{1};
        break;
        
    else
        n_var = n_var + 1;
    end 
end
deg = ones(n_var, 1).*deg;

for k=1:n_var
    name = varargin{k};
    evalin('caller', ...
        [name ' = yop.ast_control(''', name, ''',' num2str(deg(k)), ');']);
end
end