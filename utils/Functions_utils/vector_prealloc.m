classdef vector_prealloc < handle
    
    properties
        Num
        data
        id
    end
    
    methods 
        function obj = vector_prealloc(Num)
            obj.Num = Num;
            obj.data = zeros(obj.Num,1);
            obj.id = 0;
        end
        function append(obj, v)
            [r,c] = size(v);
            if c > r
                v = v';
            end
            obj.id = obj.id + length(v);
            if obj.id > obj.Num
                obj.Num = 2*obj.Num;
                obj.data = [obj.data; zeros(obj.Num,1)];
            end
            id_start = obj.id - length(v) + 1;
            id_end = obj.id;
            obj.data(id_start : id_end) = v;
        end
        function trim(obj)
            obj.data(obj.id+1:end) = [];
        end
    end
        
    
    
end