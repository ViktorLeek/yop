function ig = guess(varargin)

v=yop.ivp_var.empty(1,0);
t=v; t0=v; tf=v; x=v; z=v; u=v; p=v;
visited = [];

k=1; K=length(varargin);
while k <= K
    var_k = varargin{k};
    val_k = varargin{k+1};
    find_variables(var_k);
    [type, id] = Type(var_k);
    for r=1:length(type)
        vr = get_var(id(r));
        vr.iv = val_k(:,r);
    end
    k = k + 2;
end

if isempty(t) && (isempty(t0) || isempty(tf))
    error(yop.error.guess_missing_independent());
end

if isempty(t)
    t = yop.ivp_var(yop.independent);
    t.iv = [t0.iv; tf.iv];
end

if isscalar(t.iv)
    t.iv = [t0.iv; tf.iv];
end

if isempty(t0) && isempty(tf)
    t0 = yop.ivp_var(yop.independent0);
    tf = yop.ivp_var(yop.independentf);
    t0.iv = t.iv(1);
    tf.iv = t.iv(end);
end

for var = time_vars()
    filter(var);
end

ig = yop.ivp_sol(t.iv', horzcat(x.iv)', horzcat(z.iv, u.iv)', ...
    horzcat(p.iv)', mxargs(), variables().ids);


    function find_variables(expr)
        % Find and add all special nodes of the expression
        [tsort, N, visited] = topological_sort(expr, visited);
        
        % Find special nodes, maintain topological order
        for n=1:N
            switch class(tsort{n})
                case 'yop.ast_independent'
                    t(end+1) = yop.ivp_var(tsort{n});
                case 'yop.ast_independent_initial'
                    t0(end+1)= yop.ivp_var(tsort{n});
                case 'yop.ast_independent_final'
                    tf(end+1)= yop.ivp_var(tsort{n});
                case 'yop.ast_state'
                    x(end+1) = yop.ivp_var(tsort{n});
                case 'yop.ast_algebraic'
                    z(end+1) = yop.ivp_var(tsort{n});
                case 'yop.ast_control'
                    u(end+1) = yop.ivp_var(tsort{n});
                case 'yop.ast_parameter'
                    p(end+1) = yop.ivp_var(tsort{n});
                otherwise
                    continue
            end
        end
    end

    function v = get_var(id)
        for v = variables()
            if v.ast.m_id == id
                return;
            end
        end
        error(yop.error.failed_to_find_variable(id));
    end

    function v = variables()
        v = [t0, tf, t, x, z, u, p];
    end

    function v = time_vars()
        v = [t, x, z, u];
    end

    function mx = mxargs()
        zu = [z, u];
        mx = {t0.mx_vec(), tf.mx_vec(), t.mx_vec(), x.mx_vec, zu.mx_vec,...
            p.mx_vec()};
    end

    function filter(v)
        if isscalar(v.iv) % Scalar guess
            v.iv = interp1([t0.iv; tf.iv], [v.iv; v.iv], t.iv);
            
        elseif size(v.iv, 1) ~= length(t.iv) && size(v.iv, 1) ~= 2 
            % Guess does not confirm to scalar value, boundary values, or
            % grid.
            error(yop.error.guess_improper_grid());
            
        elseif size(v.iv, 1) ~= size(t.iv)
            % Guess at boundary values
            v.iv = interp1([t0.iv; tf.iv], v.iv, t.iv);
            
        end
    end

end