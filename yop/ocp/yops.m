function yops(varargin)
% YOPS - Declare problem variables
%   There are are seven types of variables in Yop
%       1) The independent variable, typically time
%       2) The initial value of the independent variable
%       3) The final value of the independent variable
%       4) State variables
%       5) Algebraic variables
%       6) Control inputs
%       7) Paramters
%   The initial and final value of the independent variable are parameters
%   in the optimization problems, that can either be fixed or part of the
%   OCP. Parameters are also constant variables, they can be fixed or
%   optimized for.
% 
%   EXAMPLES
%       'yops Times: t t0 tf' - declare the independent variable, its
%                               initial value and final value
%       'yops States: x1 x2 x3' - Three scalar states x1, x2 and x3
%       'yops State: x1 x2 x3' - Same as above
%       'yops State: x size: [3,1]' - A state, x, having size [3,1]
%       'yops Control: u1 u2' - Two scalar control inputs u1 and u2
%       'yops Controls: u1 u2' - Same as above with extra s for readability
%       'yops Controls: u1 u2 Control: u size: [3,1]' - Two scalar control
%                                                       inputs u1 and u2,
%                                                       and one vector
%                                                       valued control
%                                                       input u
%       'yops Ctrl: u' - Scalar control input u, less two write than above
%       'yops Parameters: p1 p2' - Two scalar parameters p1 and p2
%       'yops Params: p1 p2' - Same as above
%       'yops Param: p1 p2' - Same as above
%       'yops Param: p1 p2 Param: p3 size: [2,1]' -Two scalar parameters p1
%                                                  and p2, and one vector
%                                                  valued parameter p3
%       'yops Algebraic: z' - A scalar valued algebraic variable
%       'yops Algebraics: z' - Same as above
%       'yops Algs: z' - Same as above
%       'yops Alg: z' - Same as above
%       
%       The following line declares the independent variable (t), its
%       initial (t0) and final value (tf), two scalar states (x1, x2), a
%       vector valued control input (u) and a scalar parameter (p):
%       'yops Times: t t0 tf States: x1 x2 Control: u size: [2,1] Param: p'
%           
% 
%   INDEPENDENT VARIABLE
%   The follwing lines declares the independent variable as 't',
%   its initial value as 't0', and its final value as 'tf'
%       yops Independent: t
%       yops Independent0: t0
%       yops Independentf: tf
%   In case the independent variable is time it can be convenient to use
%   the following lines instead:
%       yops Time: t
%       yops Time0: t0
%       yops Timef: tf
%   Or you you prefer to write it on a single line
%       yops Time: t Time0: t0 Timef tf
%   Because these three variable are so common, there is a convenience call
%   to declare all three:
%      'yops Times: t t0 tf' or 'yops Independents: t t0 tf'
%   This call is parsed based on position. This means Yop will declare the 
%   first one as the independent variable, the second one as its initial 
%   value, and the third one as its final value.
%   
%   STATES
%   To declare state variables the following line can be used:
%       'yops State: x' or equivalently 'yops States: x'
%   The extra 's' does not convey a different semantic meaning here, it is
%   simply used provide slightly better readablity should you wish to
%   declare several states variables at once:
%       yops States: x1 x2 x3
%   Note that the name does not have to begin with x, it could be any valid
%   identifier name. As all states are not scalar it is possible to provide
%   a size attribute when declaring a variable:
%       yops State: x size: [3,1]
%   It is not possible to mix scalar and matrix valued states in a single
%   declaration. Instead the following can be used:
%       yops State: x size: [3,1] States: x1 x2 x3
%   This creates 4 variables, one [3,1] vector valued state 'x' and three
%   scalar states 'x1', 'x2', and 'x3'.
% 
%   ALGEBRAICS
%   Algebraic variables can be declared in the same way as states, but
%   using any of the specifiers:
%       'Algebraic:', 'Algebraics:', 'Alg:', 'Algs:'
% 
%   CONTROL INPUTS
%   Control inputs can be declared using the specifiers:
%       'Control:', 'Controls:', 'Ctrl:', 'Ctrls' 
% 
%   PARAMETERS
%   Parameters can be declared using:
%       'Parameter:', 'Parameters:', 'Param:', 'Params:'
% 
%   NOMINAL VALUES AND SCALING
%   Yop does scaling of the optimal control problem based on the user
%   input. Yop scalar variables according to:
%       v_s = (v - os)/N
%   where v_s is the scaled variable, os i the variable offset and N is the
%   nominal value of the variable. This gives that the scaling factor of
%   the shifted variable is 1/N. To be more precise N does not have to be
%   its nominal value, it could be any suitable value, but by providing
%   nominal values reasonable scaling is applied.
% 
%   CONTROL AUGMENTATION
%   One way of increasing resultion in the control input or limiting the
%   input rate, is by control its derivative and integrating that over the
%   problem horizon to obtain the control input. This means for a piecewise 
%   constant derivative that the control input is piecewise linear. How
%   many times to integrate upto the control input is specfied by adding
%   the attribute 'int:'. Example:
%       yops Control: u int: 2
%   This declares a scalar control input that is integrated from its second
%   derivative. By calling 'der(u)' the first derivaitve is obtained, and
%   by caling 'der(der(u))' the second derivative is returned.
%   
%   For vector valued control inputs the number of integrations can be
%   selected for all components by specifying a scalar value, or it can be
%   specified on a component basis as:
%       yops Control: u size: [3,1] int: [1,2,3].
%   Note that the integration value is a row vector independently of the
%   size of the control input and the enumeration of the components follows
%   that of Matlab.
t_token     = {'time:', 'independent:'};
t0_token    = {'time0:', 'independent0:'};
tf_token    = {'timef:', 'independentf:'};
ts_token    = {'times:', 'independents:'};
state_token = {'state:', 'states:'};
alg_token   = {'algebraic:', 'algebraics:', 'alg:', 'algs:'};
ctrl_token  = {'control:', 'controls:', 'ctrl:', 'ctrls:'};
param_token = {'parameter:', 'parameters:', 'param:', 'params:'};
var_flag = [t_token(:)', t0_token(:)', tf_token(:)', ts_token(:)', ...
    state_token(:)', alg_token(:)', ctrl_token(:)', param_token(:)'];
size_attr  = {'size:', 'sz:'};
w_attr     = {'nominal:', 'nom:'};
os_attr    = {'offset:', 'os:'};
int_attr   = {'int:', 'integrate:', 'aug:', 'augment:'};

decl = {};
k=1;
while k <= length(varargin)
    
    switch lower(varargin{k})
        case t_token
            step();
            time();
            
        case t0_token
            step();
            time0();
            
        case tf_token
            step();
            timef();
            
        case ts_token
            step();
            times();
            
        case state_token
            step();
            state();
            
        case alg_token
            step();
            algebraic();
            
        case ctrl_token
            step();
            control();
            
        case param_token
            step();
            parameter();
            
        otherwise
            error(yop.error.cannot_parse_yops());
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
            switch lower(varargin{k})
                case w_attr
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_attr
                    step();
                    os = str2num(varargin{k});
                    step();
                    
                case var_flag
                    break
                    
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
            switch lower(varargin{k})
                case w_attr
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_attr
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
            switch lower(varargin{k})
                case w_attr
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_attr
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
            switch lower(varargin{k})
                case w_attr
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_attr
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
            switch lower(varargin{k})
                case size_attr
                    step();
                    size = str2num(varargin{k});
                    step();
                
                case w_attr
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_attr
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
            switch lower(varargin{k})
                case size_attr
                    step();
                    size = str2num(varargin{k});
                    step();
                    
                case w_attr
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_attr
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
        vars={}; size=[1,1]; w=1; os=0; int=0;
        while k <= length(varargin)
            switch lower(varargin{k})
                case size_attr
                    step();
                    size = str2num(varargin{k});
                    step();
                    
                case w_attr
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_attr
                    step();
                    os = str2num(varargin{k});
                    step();
                    
                case var_flag
                    break;
                    
                case int_attr
                    step();
                    int = str2num(varargin{k});
                    step();
                    
                otherwise
                    vars{end+1} = varargin{k};
                    step();
            end
        end
        add_control(vars, size, w, os, int);
    end

    function parameter()
        vars={}; size=[1,1]; w=1; os=0;
        while k <= length(varargin)
            switch lower(varargin{k})
                case size_attr
                    step();
                    size = str2num(varargin{k});
                    step();
                    
                case w_attr
                    step();
                    w = str2num(varargin{k});
                    step();
                    
                case os_attr
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

    function add_control(vars, sz, w, os, int)
        if isequal(sz, [1, 1])
            add_control_scalar(vars, w, os, int);
        else
            add_control_matrix(vars, sz, w, os, int);
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

    function add_control_scalar(vars, w, os, int)
        N = length(vars);
        w   = ones(N, 1) .* w(:);
        os  = ones(N, 1) .* os(:);
        int = ones(N, 1) .* int(:);
        for n=1:N
            decl{end+1} = control_string_scalar(vars{n}, w(n), os(n), int(n));
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

    function add_control_matrix(vars, sz, w, os, int)
        decl{end+1} = control_string_matrix(vars{1}, w(:)', os(:)', sz, int);
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

function str = control_string_scalar(name, w, os, int)
str = [name, ' = yop.ast_control(''', name, ''',', num2str(w), ',', num2str(os), ',', num2str(int), ');'];
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

function str = control_string_matrix(name, w, os, sz, int)
str = [name ' = yop.control([' num2str(sz) '],''name'',''' name ''',''scaling'',[' num2str(w), '],''offset'',[' num2str(os) '], ''int'',[' num2str(int) ']);'];
end

function str = parameter_string_matrix(name, w, os, sz)
str = [name ' = yop.parameter([' num2str(sz) '],''name'',''' name ''',''scaling'',[' num2str(w), '],''offset'',[', num2str(os) ']);'];
end































