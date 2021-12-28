function yopvar(varargin)
decl = {};
t_flag     = {'-independent', '-time' , '-t'};
t0_flag    = {'-independent_initial', '-independent0', '-time0', '-t0'};
tf_flag    = {'-independent_final', '-independentf', '-timef', '-tf'};
state_flag = {'-state'};
alg_flag   = {'-algebraic', '-alg'};
ctrl_flag  = {'-control'  , '-ctrl'};
param_flag = {'-parameter', '-param'};

var_flag = [t_flag(:)', t0_flag(:)', tf_flag(:)', state_flag(:)', ...
    alg_flag(:)', ctrl_flag(:)', param_flag(:)'];

w_flag = {'-weight', '-w', '-W'};
os_flag = {'-offset', '-os', '-OS'};
deg_flag = {'-deg'};

times={}; t0={}   ; tf={}   ; states={}  ; algs={}  ; ctrls={}  ; params={}  ;
t_w=[]  ; t0_w=[] ; tf_w=[] ; state_w=[] ; alg_w=[] ; ctrl_w=[] ; param_w=[] ;
t_os=[] ; t0_os=[]; tf_os=[]; state_os=[]; alg_os=[]; ctrl_os=[]; param_os=[];
ctrl_deg = [];

k=1;
while k <= length(varargin)
    
    switch varargin{k}
        case t_flag
            step();
            time();
            
        case t0_flag
            step();
            time0();
            
        case tf_flag
            step();
            timef();
            
        case state_flag
            step();
            state();
            
        case alg_flag
            step();
            
        case ctrl_flag
            step();
            control();
            
        case param_flag
            step();
            
        otherwise
            error(yop.error.cannot_parse_yopvar());
    end

end

for r=1:length(decl)
    evalin('caller', decl{r});
end

    function step()
        k = k+1;
    end

    function time()
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    t_w = char2num(varargin{k});
                    step();
                    
                case os_flag
                    step();
                    t_os = char2num(varargin{k});
                    step();
                    
                case var_flag
                    return;
                    
                otherwise
                    times{end+1} = varargin{k};
                    step();
            end
        end
    end

    function time0()
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    t0_w = char2num(varargin{k});
                    step();
                    
                case os_flag
                    step();
                    t0_os = char2num(varargin{k});
                    step();
                    
                case var_flag
                    return;
                    
                otherwise
                    t0{end+1} = varargin{k};
                    step();
            end
        end
    end

    function timef()
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    tf_w = char2num(varargin{k});
                    step();
                    
                case os_flag
                    step();
                    tf_os = char2num(varargin{k});
                    step();
                    
                case var_flag
                    return;
                    
                otherwise
                    tf{end+1} = varargin{k};
                    step();
            end
        end
    end

    function state()
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    state_w = [state_w; char2num(varargin{k})];
                    step();
                    
                case os_flag
                    step();
                    state_os = [state_os; char2num(varargin{k})];
                    step();
                    
                case var_flag
                    return;
                    
                case deg_flag
                    
                    
                otherwise
                    states{end+1} = varargin{k};
                    step();
            end
        end
    end

    function control()
        vars={}; w=1; os=0; deg=0;
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    w = char2num(varargin{k});
                    step();
                    
                case os_flag
                    step();
                    os = char2num(varargin{k});
                    step();
                    
                case var_flag
                    break;
                    
                case deg_flag
                    step();
                    deg = char2num(varargin{k});
                    step();
                    
                otherwise
                    vars{end+1} = varargin{k};
                    step();
            end
        end
        add_control(vars, w, os, deg);
    end

    function str = add_control(vars, w, os, deg)
        N = length(vars);
        deg = ones(N, 1) .* deg;
        w   = ones(N, 1) .* w;
        os  = ones(N, 1) .* os;
        for n=1:N            
            decl{end+1} = ctrl_string(vars{n}, w(n), os(n), deg(n));
        end
    end

    function str = ctrl_string(name, w, os, deg)
        str = [name, ' = yop.ast_control(''', name, ''',', num2str(w), ',', num2str(os), ',', num2str(deg), ');'];
    end

    function num = char2num(txt)
        filt = replace(txt, ...
            {'[', ']', ', ', ' ,', ','}, ...
            { '',  '',  ' ',  ' ', ' '} ...
            );
        tmp = textscan(filt, '%f', 'delimiter', ' ');
        num = tmp{1};
    end

%%
% deg = 0;
% vars = {};
% k = 1;
% n_ctrl = 0;
% while  k <= length(varargin)
%     input = varargin{k};
%     if strcmp(input, 'deg')
%         filt = replace(varargin{k+1}, ...
%             {'[', ']', ', ', ' ,', ','}, ...
%             { '',  '',  ' ',  ' ', ' '} ...
%             );
%         tmp = textscan(filt, '%f', 'delimiter', ' ');
%         deg = tmp{1};
%         k = k+1;
%         
%     else
%         vars{end+1} = input;
%     end 
%     if input(1)=='u'
%         n_ctrl = n_ctrl + 1;
%     end
%     k = k + 1;
% end
% deg = ones(n_ctrl, 1).*deg;
% 
% cnt = 1;
% for k=1:length(vars)
%     name = vars{k};
%     if strcmp(name, 't')
%         evalin('caller', ...
%         [name ' = yop.ast_independent(''', name, ''');']);
%     
%     elseif strcmp(name, 't0')
%         evalin('caller', ...
%         [name ' = yop.ast_independent_initial(''', name, ''');']);
%     
%     elseif strcmp(name, 'tf')
%         evalin('caller', ...
%         [name ' = yop.ast_independent_final(''', name, ''');']);
%     
%     elseif name(1)=='x'
%         evalin('caller', ...
%         [name ' = yop.ast_state(''', name, ''');']);
%     
%     elseif name(1)=='z'
%         evalin('caller', ...
%         [name ' = yop.ast_algebraic(''', name, ''');']);
%     
%     elseif name(1)=='u'
%         evalin('caller', ...
%         [name ' = yop.ast_control(''', name, ''',' num2str(deg(cnt)), ');']);
%         cnt = cnt + 1;
%     
%     elseif name(1)=='p'
%         evalin('caller', ...
%         [name ' = yop.ast_parameter(''', name, ''');']);
%     
%     
%     else 
%         error(yop.error.yopvar_not_recognized(name));
%     end
% end
end