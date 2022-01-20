%% Kvar att göra
% -[ ] Persistant-variabel som väntar på att options ska sparas
% -[ ] evalin

function varargout = yoptions(varargin)
ival_flags = {'Interval:', 'interval:', 'Intervals:', 'intervals:' 'Ival:', 'ival:', 'N:', 'n:'};
state_flags = {'State:', 'state:', 'X:', 'x:'};
control_flags = {'Control:', 'control:', 'Ctrl:', 'ctrl:', 'U:', 'u:'};
solver_flags = {'Solver:', 'solver:', 'S:', 's:'};
decl_flags = {'->', '>>', '=:', '#', '.', '@'};
decl_default_flags = {'!'};
name_flags = {decl_flags{:}};
points_flags = {'legendre', 'radau'};

flags = {ival_flags{:}, state_flags{:}, control_flags{:}, ...
    solver_flags{:}, decl_flags{:}, name_flags{:}, decl_default_flags{:}};


yopts = struct;
yopts.solver = '';
yopts.intervals = [];
yopts.state_points = {};
yopts.state_deg = [];
yopts.control_points = {};
yopts.control_deg = [];

declvar = '';
k=1;
while k <= length(varargin)
    
    tmp = varargin{k};
    cont = false;
    for n=1:length(name_flags)
        nf = name_flags{n};
        if length(tmp)>=length(nf) && ~strcmp(tmp, nf) && strcmp(tmp(1:length(nf)), nf)
            declvar = tmp(length(nf)+1:end);
            step();
            cont = true;
        end
    end
    if cont
        continue
    end
    
    switch varargin{k}
        
        case decl_flags
            step();
            declare();
            
        case decl_default_flags
            step();
            declare_default();
        
        case solver_flags
            step();
            solver();
        
        case ival_flags
            step();
            interval();
            
        case state_flags
            step();
            state();
            
        case control_flags
            step();
            control();
            
        otherwise
            error(yop.error.cannot_parse_yops());
    end
    
end

yopts.declvar = declvar;
varargout{1} = yopts;
% for r=1:length(declopts)
%     evalin('caller', declopts{r});
% end

    function step()
        k=k+1;
    end

    function declare()
        declvar = varargin{k};
        step();
    end

    function declare_default()
        declvar = 'yopts';
    end

    function solver()
        yopts.solver = lower(varargin{k});
        step();
    end

    function interval()
        while k <= length(varargin)
            
            tmp = varargin{k};
            for nn=1:length(name_flags)
                nf = name_flags{nn};
                if length(tmp)>=length(nf) && ~strcmp(tmp, nf) && strcmp(tmp(1:length(nf)), nf)
                    return
                end
            end
            
            switch varargin{k}
                case flags
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
            
            tmp = varargin{k};
            for nn=1:length(name_flags)
                nf = name_flags{nn};
                if length(tmp)>=length(nf) && ~strcmp(tmp, nf) && strcmp(tmp(1:length(nf)), nf)
                    return
                end
            end
            
            switch varargin{k}
                case flags
                    break
                    
                case points_flags
                    yopts.state_points{end+1} = lower(varargin{k});
                    step();
                    
                otherwise
                    yopts.state_deg = [yopts.state_deg, ...
                        str2num(varargin{k})];
                    step();
                    
            end
        end
    end

    function control()
        while k <= length(varargin)
            
            tmp = varargin{k};
            for nn=1:length(name_flags)
                nf = name_flags{nn};
                if length(tmp)>=length(nf) && ~strcmp(tmp, nf) && strcmp(tmp(1:length(nf)), nf)
                    return
                end
            end
            
            switch varargin{k}
                case flags
                    break
                    
                case points_flags
                    yopts.control_points{end+1} = lower(varargin{k});
                    step();
                    
                otherwise
                    yopts.control_deg = [yopts.control_deg, ...
                        str2num(varargin{k})];
                    step();
                    
            end
        end
    end

end