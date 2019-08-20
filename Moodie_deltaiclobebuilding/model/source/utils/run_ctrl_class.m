classdef run_ctrl_class
    %RUN_CTRL 
    % this is the object that handles run control for the model
    % basically, while the object holds a status of "true", the run continues
    %
    
    properties
        nYears = NaN;
        nAvuls = NaN;
        counter = 0;
        incr = 1;
        st = true; % status value
        str;
        limit;
        nDays;
        day_incr = false
        avul_incr = false
        Qw_num;
        Qw_func;
    end
    
    methods
        
        function obj = run_ctrl_class(str, n)
            %RUN_CTRL Construct an instance of this class
            %   Detailed explanation goes here
            
            if strcmp(str, 'nYears')
                obj.str = str;
                obj.nYears = n;
                obj.nDays = obj.nYears * 365;
                obj.limit = obj.nDays;
                obj.day_incr = true;
                obj.Qw_func = (@(Qw_num) Qw_num + 1);
                
            elseif strcmp(str, 'nAvuls')
                obj.str = str;
                obj.nAvuls = n;
                obj.limit = obj.nAvuls;
                obj.avul_incr = true;
                obj.Qw_func = (@(Qw_num) Qw_num);
                
            else
                error('Bad run_ctrl.str supplied')
                
            end
            
        end
        
        function st = status(obj)
            %STATUS Query the state of the run_ctrl object
            %   Returns boolean on whether counter < limit
            
            st = obj.st;
            
        end
        
        function obj = increment(obj)
            %INCREMENT Increment the counter
            obj.counter = obj.counter + obj.incr;
            obj.st = obj.counter < obj.limit;
        end
        
        function obj = kill(obj)
            %kill sets the status to false to end processing
            
            obj.st = false;
            
        end
        
    end
end

