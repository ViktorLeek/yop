function yopvar(varargin)
decl = {};
t_flag     = {'-independent', '-time' , '-t', 'time:'};
t0_flag    = {'-independent_initial', '-independent0', '-time0', '-t0', 'time0:', 'time_0:'};
tf_flag    = {'-independent_final', '-independentf', '-timef', '-tf', 'timef:', 'time_f:'};
ts_flag    = {'times:'};
state_flag = {'-state', 'states:', 'state:'};
alg_flag   = {'-algebraic', '-alg', 'algs:', 'algebraics:'};
ctrl_flag  = {'-control'  , '-ctrl', 'ctrl:', 'ctrls:', 'controls:'};
param_flag = {'-parameter', '-param', 'param:', 'params:', 'parameters:'};
var_flag = [t_flag(:)', t0_flag(:)', tf_flag(:)', ts_flag(:), ...
    state_flag(:)', alg_flag(:)', ctrl_flag(:)', param_flag(:)'];

size_flag  = {'size:', 'sz:'};
w_flag = {'-weight', '-w', '-W', 'weight:', 'w:', 'W:', 'weights:', 'scaling:'};
os_flag = {'-offset', '-os', '-OS', 'offset:', 'os:', 'OS:', 'offsets:'};
deg_flag = {'-deg', 'deg:', 'degree:'};

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
            
        case ts_flag
            step();
            times();
            
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
        k=k+1;
    end

    function time()
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
        add_time(vars, w, os);
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

    function times()
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
        add_times(vars, w, os);
    end

    function state()
        vars={}; size=[1,1]; w=1; os=0;
        while k <= length(varargin)
            switch varargin{k}
                case size_flag
                    step();
                    size = str2num(varargin{k});
                    step();
                
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
        add_state(vars, size, w, os);
    end

    function algebraic()
        vars={};  size=[1,1]; w=1; os=0;
        while k <= length(varargin)
            switch varargin{k}
                case size_flag
                    step();
                    size = str2num(varargin{k});
                    step();
                    
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
        add_algebraic(vars, size, w, os);
    end

    function control()
        vars={}; size=[1,1]; w=1; os=0; deg=0;
        while k <= length(varargin)
            switch varargin{k}
                case size_flag
                    step();
                    size = str2num(varargin{k});
                    step();
                    
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
        add_control(vars, size, w, os, deg);
    end

    function parameter()
        vars={}; size=[1,1]; w=1; os=0;
        while k <= length(varargin)
            switch varargin{k}
                case size_flag
                    step();
                    size = str2num(varargin{k});
                    step();
                    
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
        add_parameter(vars, size, w, os);
    end

    function add_time(vars, w, os)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        for n=1:N
            decl{end+1} = time_string(vars{n}, w(n), os(n));
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

    function add_times(vars, w, os)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        
        if N >= 1
            decl{end+1} = time_string(vars{1}, w(1), os(1));
        end
        
        if N >= 2
            decl{end+1} = time0_string(vars{2}, w(2), os(2));
        end
        
        if N >= 3
            decl{end+1} = timef_string(vars{3}, w(3), os(3));
        end
        
    end

    function add_state(vars, sz, w, os)
        if isequal(sz, [1, 1])
            add_state_scalar(vars, w, os);
        else
            add_state_matrix(vars, sz, w, os);
        end
    end

    function add_algebraic(vars, sz, w, os)
        if isequal(sz, [1, 1])
            add_algebraic_scalar(vars, w, os);
        else
            add_algebraic_matrix(vars, sz, w, os);
        end
    end

    function add_control(vars, sz, w, os, deg)
        if isequal(sz, [1, 1])
            add_control_scalar(vars, w, os, deg);
        else
            add_control_matrix(vars, sz, w, os, deg);
        end
    end

    function add_parameter(vars, sz, w, os)
        if isequal(sz, [1, 1])
            add_parameter_scalar(vars, w, os);
        else
            add_parameter_matrix(vars, sz, w, os);
        end
    end

    function add_state_scalar(vars, w, os)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        for n=1:N
            decl{end+1} = state_string_scalar(vars{n}, w(n), os(n));
        end
    end

    function add_algebraic_scalar(vars, w, os)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        for n=1:N
            decl{end+1} = algebraic_string_scalar(vars{n}, w(n), os(n));
        end
    end

    function add_control_scalar(vars, w, os, deg)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        deg = ones(N, 1) .* deg(:);
        for n=1:N
            decl{end+1} = control_string_scalar(vars{n}, w(n), os(n), deg(n));
        end
    end

    function add_parameter_scalar(vars, w, os)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        for n=1:N
            decl{end+1} = parameter_string_scalar(vars{n}, w(n), os(n));
        end
    end

    function add_state_matrix(vars, sz, w, os)
        decl{end+1} = state_string_matrix(vars{1}, w(:)', os(:)', sz);
    end


    function add_algebraic_matrix(vars, sz, w, os)
        decl{end+1} = algebraic_string_matrix(vars{1}, w(:)', os(:)', sz);
    end

    function add_control_matrix(vars, sz, w, os, deg)
        decl{end+1} = control_string_matrix(vars{1}, w(:)', os(:)', sz, deg);
    end

end

function str = time_string(name, w, os)
str = [name, ' = yop.ast_independent(''', name, ''',', num2str(w), ',', num2str(os), ');'];
end

function str = time0_string(name, w, os)
str = [name, ' = yop.ast_independent_initial(''', name, ''',', num2str(w), ',', num2str(os), ');'];
end

function str = timef_string(name, w, os)
str = [name, ' = yop.ast_independent_final(''', name, ''',', num2str(w), ',', num2str(os), ');'];
end

function str = state_string_scalar(name, w, os)
str = [name, ' = yop.ast_state(''', name, ''',', num2str(w), ',', num2str(os) ');'];
end

function str = algebraic_string_scalar(name, w, os)
str = [name, ' = yop.ast_algebraic(''', name, ''',', num2str(w), ',', num2str(os), ');'];
end

function str = control_string_scalar(name, w, os, deg)
str = [name, ' = yop.ast_control(''', name, ''',', num2str(w), ',', num2str(os), ',', num2str(deg), ');'];
end

function str = parameter_string_scalar(name, w, os)
str = [name, ' = yop.ast_parameter(''', name, ''',', num2str(w), ',', num2str(os), ');'];
end

function str = state_string_matrix(name, w, os, sz)
str = [name ' = yop.state([' num2str(sz) '],''name'',''' name ''',''scaling'',[' num2str(w), '],''offset'',[', num2str(os) ']);'];
end

function str = algebraic_string_matrix(name, w, os, sz)
str = [name ' = yop.algebraic([' num2str(sz) '],''name'',''' name ''',''scaling'',[' num2str(w), '],''offset'',[', num2str(os) ']);'];
end

function str = control_string_matrix(name, w, os, sz, deg)
str = [name ' = yop.control([' num2str(sz) '],''name'',''' name ''',''scaling'',[' num2str(w), '],''offset'',[' num2str(os) '], ''deg'',[' num2str(deg) ']);'];
end

function str = parameter_string_matrix(name, w, os, sz)
str = [name ' = yop.parameter([' num2str(sz) '],''name'',''' name ''',''scaling'',[' num2str(w), '],''offset'',[', num2str(os) ']);'];
end































