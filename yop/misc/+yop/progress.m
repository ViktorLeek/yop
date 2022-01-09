classdef progress
    methods (Static)
        
        function ocp_parsing()
            fprintf("[Yop]: Parsing optimal control problem.\n");
        end
        
        function ocp_parsed()
            fprintf("[Yop]: Optimal control problem parsed.\n");
        end
        
        function ocp_solved(solved)
            if solved
                fprintf("[Yop]: Yet another yoptimal control problem solved.\n");
            else
                fprintf("[Yop]: Optimal control problem failed to converge.\n");
            end
        end
        
        function msg = nlp_building_msg()
            msg = "[Yop]: Building NLP: [________]";
        end
        
        function msg = nlp_initialized_msg()
            msg = "[Yop]: Building NLP: [#_______]";
        end
        
        function msg = nlp_special_nodes_done_msg()
            msg = "[Yop]: Building NLP: [##______]";
        end
        
        function msg = nlp_dynamics_done_msg()
            msg = "[Yop]: Building NLP: [###_____]";
        end
        
        function msg = nlp_pointcons_done_msg()
            msg = "[Yop]: Building NLP: [####____]";
        end
        
        function msg = nlp_pathcons_done_msg()
            msg = "[Yop]: Building NLP: [#####___]";
        end
        
        function msg = nlp_box_done_msg()
            msg = "[Yop]: Building NLP: [######__]";
        end
        
        function msg = nlp_guess_done_msg()
            msg = "[Yop]: Building NLP: [#######_]";
        end
        
        function msg = nlp_completed_msg()
            msg = "[Yop]: Building NLP: [########]\n";
        end
        
        function nlp_building()
            fprintf(yop.progress.nlp_building_msg());
        end
        
        function nlp_initialized()
            yop.progress.erase(yop.progress.nlp_building_msg());
            fprintf(yop.progress.nlp_initialized_msg());
        end
        
        function nlp_special_nodes_done()
            yop.progress.erase(yop.progress.nlp_initialized_msg());
            fprintf(yop.progress.nlp_special_nodes_done_msg());
        end
        
        function nlp_dynamics_done()
            yop.progress.erase(yop.progress.nlp_special_nodes_done_msg());
            fprintf(yop.progress.nlp_dynamics_done_msg());
        end
        
        function nlp_pointcons_done()
            yop.progress.erase(yop.progress.nlp_dynamics_done_msg());
            fprintf(yop.progress.nlp_pointcons_done_msg());
        end
        
        function nlp_pathcons_done()
            yop.progress.erase(yop.progress.nlp_pointcons_done_msg());
            fprintf(yop.progress.nlp_pathcons_done_msg());
        end
        
        function nlp_box_done()
            yop.progress.erase(yop.progress.nlp_pathcons_done_msg());
            fprintf(yop.progress.nlp_box_done_msg());
        end
        
        function nlp_guess_done()
            yop.progress.erase(yop.progress.nlp_box_done_msg());
            fprintf(yop.progress.nlp_guess_done_msg());
        end
        
        function nlp_completed()
            yop.progress.erase(yop.progress.nlp_guess_done_msg());
            fprintf(yop.progress.nlp_completed_msg());
        end
        
        function erase(msg)
            fprintf(repmat('\b', 1, strlength(msg)));
        end
    end
end