function varargout = yoptions(varargin)
% YOPTIONS - Create options to Yop
% This is a very flexible function for creating options which also means
% that it comes with side effects which must be considered. 
% 
% Everytime you read the options that has been entered they are consumed. 
% You can read out default options without consuming the input options 
% by writing: 
%   yoptions('default') 
% 
% The options are returned by calling yoptions without arguments
%   opts = yoptions
% or at the end of the options entering any of the commands:
%   '-> optsname' 
%   '>> optsname', 
%   '@optsname' (no space!) 
% they are the save to a variable name 'optsname' which is could be any
% valid identifier name. The can also be returned to a by default name
% variable yopts by entering any of the commands:
%   '.'
%   '!' 
% 
% YOPTIONS - Example use:
%   yoptions Ivals: 100 State: 5 L Control: 2 R 
%   yoptions -ipopt -continuity max_iter: 5000 -> yopts
%   
%   The first line enters the options:
%       Ivals: 100 - 100 control intervals
%       State: 5 L - state polynomial degree 5, Legendre collocation points
%       Control: 2 R - control poly degree 2, Radau collocation points
% 
%       Since there is no command for returning the options are given they
%       are remembered.
% 
%   The second line enters the options
%       -ipopt - which selects IPOPT as solver. 
%       -continuity  - which means control input is continuous (relevant
%                      only when the control degree > 0)
%       max_iter: 5000 - Sets the maximum number of ipopt solver
%                         iterations to 5000. This does not have to be
%                         preceded by the "-ipopt" option. Although it is
%                         ignored if another solver is selected.
%       -> yopts - which means to return the options to a variable named
%                  yopts. The name could be any valid variable name.

persistent yopts
if isempty(yopts)
    reset();
end

if nargin == 0
    insert_defaults();
    varargout{1} = yopts;
    reset();
    return
    
elseif nargin == 1 && strcmp(varargin{1}, 'default')
    varargout{1} = get_defaults();
    return
end

% Top level tokens
ival_tokens = { ...
    'interval:', ...
    'intervals:', ...
    'ival:', ...
    'ivals:', ...
    'n:' ...
    };
state_tokens = { ...
    'state:', ...
    'x:' ...
    };
control_tokens = { ...
    'control:', ...
    'ctrl:', ...
    'u:' ...
    };
solver_tokens = { ...
    'solver:', ...
    's:' ...
    };
decl_tokens = { ...
    '->', ...
    '>>', ...
    '@' ...
    };
decl_default_tokens = { ...
    '!', ...
    '.', ...
    '-decl', ...
    '-declare', ...
    '-make', ...
    '-save' ...
    };
name_tokens = { ...
    decl_tokens{:} ...
    };
reset_tokens = { ...
    '~', ...
    '-reset', ...
    '-clc', ...
    '-clear' ...
    };
continuity_tokens = { ...
    '-continuity', ...
    '-continuous_control', ...
    '-control_continuity', ...
    '-cc', ...
    'continuity:', ...
    '-link' ...
    };
ipopt_tokens = { ...
    '-ipopt' ...
    };
max_iter_tokens = { ...
    'max_iter:', ...
    'maxiter:', ...
    'max_it:', ...
    'maxit:' ...
    };
acc_iter_tokens = { ...
    'acceptable_iter:', ...
    'acc_iter:', ...
    'aiter:', ...
    'acciter:' ...
    };
acctol_tokens = { ...
    'atol:', ...
    'acctol:', ...
    'acceptable_tol:' ...
    };
opttol_tokens = { ...
    'tol:', ...
    'opt_tol:', ...
    'optimality_tol' ...
    };
constr_viol_tokens = { ...
    'constr_viol_tol:', ...
    'constr_tol:', ...
    'cvtol:', ...
    'ctol:' ...
    };
print_level_tokens = { ...
    'print_level:', ...
    'print:', ...
    'plevel:' ...
    };
linear_solver_tokens = { ...
    'linear_solver:', ...
    'lin_sol:', ...
    'ls:' ...
    };
ipopt_options_tokens = { ...
    max_iter_tokens{:}, ...
    acc_iter_tokens{:}, ...
    acctol_tokens{:}, ...
    constr_viol_tokens{:}, ...
    print_level_tokens{:}, ...
    linear_solver_tokens{:}, ...
    opttol_tokens{:} ...
    };

top_lvl_tokens = {ival_tokens{:}, state_tokens{:}, control_tokens{:}, ...
    solver_tokens{:}, decl_tokens{:}, decl_default_tokens{:}, name_tokens{:},...
    reset_tokens{:}, continuity_tokens{:}, ipopt_options_tokens{:}, ...
    ipopt_tokens{:}};

% State specific tokens
legendre_tokens = { ...
    'legendre', ...
    'l', ...
    '-l' ...
    };
radau_tokens = { ...
    'radau', ...
    'r', ...
    '-r' ...
    };
points_tokens = { ...
    legendre_tokens{:}, ...
    radau_tokens{:} ...
    };

declname = '';
declarations = {};

k=1;
name_ptr = yop.pointer();
while k <= length(varargin)
    
    if toplevel_declaration(name_ptr)
        declname = name_ptr.value;
        step();
        add_declaration();
        continue
    end
    
    switch lower(varargin{k})
        
        case decl_tokens
            step();
            declare();
            
        case decl_default_tokens
            step();
            declare_default();
            
        case continuity_tokens
            step();
            continuity();
            
        case reset_tokens
            step();
            reset();
        
        case solver_tokens
            step();
            solver();
        
        case ival_tokens
            step();
            interval();
            
        case state_tokens
            step();
            state();
            
        case control_tokens
            step();
            control();
            
        case ipopt_tokens
            step();
            set_solver('ipopt');
            
        case max_iter_tokens
            step();
            max_iter();
            
        case acc_iter_tokens
            step();
            acc_iter();
            
        case acctol_tokens
            step();
            acc_tol();
            
        case opttol_tokens
            step();
            opt_tol();
            
        case constr_viol_tokens
            step();
            constr_viol();
            
        case print_level_tokens
            step();
            print_level();

        case linear_solver_tokens
            step();
            linear_solver();
            
        otherwise
            error(yop.error.cannot_parse_yops());
    end
    
end

for r=1:length(declarations)
    evalin('caller', declarations{r});
end

    function step()
        k=k+1;
    end

    function reset()
        yopts = struct;
        yopts.solver         = '';
        yopts.intervals      = [];
        yopts.state_points   = {};
        yopts.state_degree   = [];
        yopts.control_points = {};
        yopts.control_degree = [];
        yopts.continuity     = [];
        yopts.ipopts         = struct;
    end
    
    function insert_defaults()
        if isempty(yopts.solver)
            yopts.solver = yop.defaults.solver;
        end
        if isempty(yopts.intervals)
            yopts.intervals = yop.defaults.control_invervals;
        end
        if isempty(yopts.state_points)
            yopts.state_points = yop.defaults.state_points;
        end
        if isempty(yopts.state_degree)
            yopts.state_degree = yop.defaults.state_degree;
        end
        if isempty(yopts.control_points)
            yopts.control_points = yop.defaults.control_points;
        end
        if isempty(yopts.control_degree)
            yopts.control_degree = yop.defaults.control_degree;
        end
        if isempty(yopts.continuity)
            yopts.continuity = yop.defaults.continuity;
        end
    end

    function opts = get_defaults()
        opts.solver = yop.defaults.solver;
        opts.intervals = yop.defaults.control_invervals;
        opts.state_points = yop.defaults.state_points;
        opts.state_degree = yop.defaults.state_degree;
        opts.control_points = yop.defaults.control_points;
        opts.control_degree = yop.defaults.control_degree;
        opts.continuity = yop.defaults.continuity;
        opts.ipopts = struct;
    end

    function declare()
        declname = varargin{k};
        step();
        add_declaration()
    end

    function declare_default()
        declname = 'yopts';
        add_declaration()
    end

    function continuity()
        yopts.continuity = true;
    end

    function solver()
        set_solver(lower(varargin{k}));
        step();
    end

    function set_solver(name)
        yopts.solver = name;
    end

    function max_iter()
        yopts.ipopts.max_iter = str2num(varargin{k});
        step();
    end

    function acc_iter()
        yopts.ipopts.acceptable_iter = str2num(varargin{k});
        step();
    end

    function acc_tol()
        yopts.ipopts.acceptable_tol = str2num(varargin{k});
        step();
    end

    function opt_tol()
        yopts.ipopts.tol = str2num(varargin{k});
        step();
    end

    function constr_viol()
        yopts.ipopts.constr_viol_tol = str2num(varargin{k});
        step();
    end

    function print_level()
        yopts.ipopts.print_level = str2num(varargin{k});
        step();
    end

    function linear_solver()
        yopts.ipopts.linear_solver = varargin{k};
        step();
    end

    function interval()
        while k <= length(varargin)
            
            if toplevel_declaration()
                return
            end
            
            switch lower(varargin{k})
                case top_lvl_tokens
                    break
                    
                otherwise
                    yopts.intervals = [yopts.intervals, ...
                        str2num(varargin{k})];
                    step();
                    
            end
        end
    end

    function state()
        while k <= length(varargin)
            
            if toplevel_declaration()
                return
            end
            
            switch lower(varargin{k})
                case top_lvl_tokens
                    break
                    
                case points_tokens
                    switch lower(varargin{k})
                        case legendre_tokens
                            pnts = 'legendre';
                        case radau_tokens
                            pnts = 'radau';
                    end
                    yopts.state_points{end+1} = pnts;
                    step();
                    
                otherwise
                    yopts.state_degree = [yopts.state_degree, ...
                        str2num(varargin{k})];
                    step();
                    
            end
        end
    end

    function control()
        while k <= length(varargin)
            
            if toplevel_declaration()
                return
            end
            
            switch lower(varargin{k})
                case top_lvl_tokens
                    break
                    
                case points_tokens
                    switch lower(varargin{k})
                        case legendre_tokens
                            pnts = 'legendre';
                        case radau_tokens
                            pnts = 'radau';
                    end
                    yopts.control_points{end+1} = pnts;
                    step();
                    
                otherwise
                    yopts.control_degree = [yopts.control_degree, ...
                        str2num(varargin{k})];
                    step();
                    
            end
        end
    end

    function bool = toplevel_declaration(pointer)
        bool = false;
        cur = varargin{k};
        for nn=1:length(name_tokens)
            nf = name_tokens{nn};
            if length(cur)>=length(nf) && ~strcmp(cur, nf) && strcmp(cur(1:length(nf)), nf)
                bool = true;
                if nargin == 1
                    name = cur(length(nf)+1:end);
                    pointer.value = name;
                end
            end
        end
    end

    function add_declaration()
        insert_defaults();
        str = [declname '=struct('];
        fields = fieldnames(yopts);
        N  = numel(fields);
        for n=1:N
            field = fields{n};
            value = yopts.(field);
            str = [str, q(field), ','];
            if isempty(value)
                str = [str, '[]'];
            elseif isnumeric(value) || islogical(value)
                str = [str, b(num2str(value))];
            elseif iscellstr(value)
                str = [str, cellstr2str(value)];
            elseif ischar(value) || isstring(value)
                str = [str, q(value)];
            elseif isstruct(value) && strcmp(field, 'ipopts')
                str = [str, ipopt2str(value)];
            else
                error(yop.error.unexpected_error());
            end
            str = yop.IF(n==N, str, [str, ',']);
        end
        str = [str, ');'];
        declarations{end+1} = str;
        reset();
    end

    function str = ipopt2str(ipopt)
        str = 'struct(';
        fields = fieldnames(ipopt);
        N  = numel(fields);
        for n=1:N
            field = fields{n};
            value = ipopt.(field);
            if isstring(value) || ischar(value)
                str = [str, q(field), ',' q(value)];
            else
                str = [str, q(field), ',' num2str(value)];
            end
            str = yop.IF(n==N, str, [str, ',']);
        end
        str = [str, ')'];
    end

    function qted = q(f)
        qted = ['''', f, ''''];
    end

    function bted = b(numstr)
        bted = ['[',numstr,']'];
    end

    function c = cellstr2str(cs)
        if isempty(cs)
            c = '{{}}';
            return;
        end
        c = '{{';
        for n=1:length(cs)-1
            c = [c, q(cs{n}), ','];
        end
        c = [c, q(cs{end}), '}}'];
    end

end