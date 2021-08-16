classdef ast_sum < yop.ast_expression
    properties
        expr
        opt1
        opt2
        opt3
        nargs
    end
    methods
        function obj = ast_sum(expr, opt1, opt2, opt3)
            % ast_sum
            
            % This only works for cases where the optional arguments are
            % not ast_node
            obj.expr = expr;
            obj.nargs = nargin;
            switch nargin
                case 1
                    obj.dim = size(sum(ones(size(expr))));
                    
                case 2
                    obj.opt1 = opt1;
                    obj.dim = size(sum(ones(size(expr)), opt1));
                    
                case 3
                    obj.opt1 = opt1;
                    obj.opt2 = opt2;
                    obj.dim = size(sum(ones(size(expr)), opt1, opt2));
                    
                case 4                    
                    obj.opt1 = opt1;
                    obj.opt2 = opt2;
                    obj.opt3 = opt3;
                    obj.dim = size(sum(ones(size(expr)), opt1, opt2, opt3));
            end
        end
        
        function value = evaluate(obj)
            switch obj.nargs
                case 1
                    value = sum(evaluate(obj.expr));
                    
                case 2
                    value = sum(evaluate(obj.expr), obj.opt1);
                    
                case 3
                    value = sum(evaluate(obj.expr), obj.opt1, obj.opt2);
                    
                case 4
                    value = sum(...
                        evaluate(obj.expr), ...
                        obj.opt1, ...
                        obj.opt2, ...
                        obj.opt3 ...
                        );
                    
            end
        end
        
        function draw(obj)
            
            switch obj.nargs
                case 1
                    fprintf('sum(expr)\n');
                    last_child(obj);
                    draw(obj.expr);
                    end_child(obj);
                    
                case 2
                    fprintf('sum(expr, opt1)\n');
                    
                    begin_child(obj);
                    draw(obj.expr);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.opt1);
                    end_child(obj);
                    
                case 3
                    fprintf('sum(expr, opt1, opt2)\n');
                    
                    begin_child(obj);
                    draw(obj.expr);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.opt1);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.opt2);
                    end_child(obj);
                    
                case 4
                    fprintf('sum(expr, opt1, opt2, opt3)\n');
                    
                    begin_child(obj);
                    draw(obj.expr);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.opt1);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.opt2);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.opt3);
                    end_child(obj);
            end
            
        end
    end
end