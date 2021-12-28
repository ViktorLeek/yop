function yopvar(varargin)
decl = {};
t_flag     = {'-independent', '-time' , '-t', 'time:'};
t0_flag    = {'-independent_initial', '-independent0', '-time0', '-t0', 'time0:'};
tf_flag    = {'-independent_final', '-independentf', '-timef', '-tf', 'timef:'};
state_flag = {'-state', 'states:'};
alg_flag   = {'-algebraic', '-alg', 'algs:', 'algebraics:'};
ctrl_flag  = {'-control'  , '-ctrl', 'ctrls:', 'controls:'};
param_flag = {'-parameter', '-param', 'params:', 'parameters:'};

var_flag = [t_flag(:)', t0_flag(:)', tf_flag(:)', state_flag(:)', ...
    alg_flag(:)', ctrl_flag(:)', param_flag(:)'];

w_flag = {'-weight', '-w', '-W', 'weight:', 'w:', 'W:', 'weights:', 'scaling:'};
os_flag = {'-offset', '-os', '-OS', 'offset:', 'os:', 'OS:', 'offsets:'};
deg_flag = {'-deg', 'deg:'};

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
            algebraic();
            
        case ctrl_flag
            step();
            control();
            
        case param_flag
            step();
            parameter();
            
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
        vars={};
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    error(yop.error.independent_scaled());
                    
                case os_flag
                    error(yop.error.independent_offset());
                    
                case var_flag
                    break;
                    
                otherwise
                    vars{end+1} = varargin{k};
                    step();
            end
        end
        add_time(vars);
    end

    function time0()
        vars={}; w=1; os=0;
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_flag
                    step();
                    os = str2num(varargin{k});
                    step();
                    
                case var_flag
                    break;
                    
                otherwise
                    vars{end+1} = varargin{k};
                    step();
            end
        end
        add_time0(vars, w, os);
    end

    function timef()
        vars={}; w=1; os=0;
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_flag
                    step();
                    os = str2num(varargin{k});
                    step();
                    
                case var_flag
                    break;
                    
                otherwise
                    vars{end+1} = varargin{k};
                    step();
            end
        end
        add_timef(vars, w, os);
    end

    function state()
        vars={}; w=1; os=0;
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_flag
                    step();
                    os = str2num(varargin{k});
                    step();
                    
                case var_flag
                    break;
                    
                otherwise
                    vars{end+1} = varargin{k};
                    step();
            end
        end
        add_state(vars, w, os);
    end

    function algebraic()
        vars={}; w=1; os=0;
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_flag
                    step();
                    os = str2num(varargin{k});
                    step();
                    
                case var_flag
                    break;
                    
                otherwise
                    vars{end+1} = varargin{k};
                    step();
            end
        end
        add_algebraic(vars, w, os);
    end

    function control()
        vars={}; w=1; os=0; deg=0;
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_flag
                    step();
                    os = str2num(varargin{k});
                    step();
                    
                case var_flag
                    break;
                    
                case deg_flag
                    step();
                    deg = str2num(varargin{k});
                    step();
                    
                otherwise
                    vars{end+1} = varargin{k};
                    step();
            end
        end
        add_control(vars, w, os, deg);
    end

    function parameter()
        vars={}; w=1; os=0;
        while k <= length(varargin)
            switch varargin{k}
                case w_flag
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_flag
                    step();
                    os = str2num(varargin{k});
                    step();
                    
                case var_flag
                    break;
                    
                otherwise
                    vars{end+1} = varargin{k};
                    step();
            end
        end
        add_parameter(vars, w, os);
    end

    function add_time(vars)
        N = length(vars);
        for n=1:N
            decl{end+1} = time_string(vars{n});
        end
    end

    function add_time0(vars, w, os)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        for n=1:N
            decl{end+1} = time0_string(vars{n}, w(n), os(n));
        end
    end

    function add_timef(vars, w, os)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        for n=1:N
            decl{end+1} = timef_string(vars{n}, w(n), os(n));
        end
    end

    function add_state(vars, w, os)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        for n=1:N
            decl{end+1} = state_string(vars{n}, w(n), os(n));
        end
    end

    function add_algebraic(vars, w, os)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        for n=1:N
            decl{end+1} = algebraic_string(vars{n}, w(n), os(n));
        end
    end

    function add_control(vars, w, os, deg)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        deg = ones(N, 1) .* deg(:);
        for n=1:N            
            decl{end+1} = ctrl_string(vars{n}, w(n), os(n), deg(n));
        end
    end

    function add_parameter(vars, w, os)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        for n=1:N
            decl{end+1} = parameter_string(vars{n}, w(n), os(n));
        end
    end

    function str = time_string(name)
        str = [name, ' = yop.ast_independent(''', name, ''');'];
    end

    function str = time0_string(name, w, os)
        str = [name, ' = yop.ast_independent_initial(''', name, ''',', num2str(w), ',', num2str(os), ');'];
    end

    function str = timef_string(name, w, os)
        str = [name, ' = yop.ast_independent_final(''', name, ''',', num2str(w), ',', num2str(os), ');'];
    end

    function str = state_string(name, w, os)
        str = [name, ' = yop.ast_state(''', name, ''',', num2str(w), ',', num2str(os), ');'];
    end

    function str = algebraic_string(name, w, os)
        str = [name, ' = yop.ast_algebraic(''', name, ''',', num2str(w), ',', num2str(os), ');'];
    end

    function str = ctrl_string(name, w, os, deg)
        str = [name, ' = yop.ast_control(''', name, ''',', num2str(w), ',', num2str(os), ',', num2str(deg), ');'];
    end

    function str = parameter_string(name, w, os)
        str = [name, ' = yop.ast_parameter(''', name, ''',', num2str(w), ',', num2str(os), ');'];
    end

%     function num = char2num(txt)
%         filt = replace(txt, ...
%             {'[', ']', ', ', ' ,', ','}, ...
%             { '',  '',  ' ',  ' ', ' '} ...
%             );
%         tmp = textscan(filt, '%f', 'delimiter', ' ');
%         num = tmp{1};
%     end

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