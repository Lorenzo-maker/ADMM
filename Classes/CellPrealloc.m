classdef CellPrealloc < handle

    properties
        Data
        index
        sz
    end


    methods

        function obj = CellPrealloc(sz)
            
            % init cell
            obj.Data = cell(sz, 1);
            obj.sz = sz;
            obj.index = 0;
        end

        function append(obj, appendData)
            
            [r, c] = size(appendData);
            if c>r
                appendData = appendData';
            end
            
            obj.index = obj.index+1;
            if obj.index > obj.sz
                obj.realloc();
            end
            obj.Data{obj.index} = appendData;

        end
        
        function realloc(obj)
            obj.Data = [obj.Data; cell(obj.sz, 1)];
        end

        function trim(obj)
            obj.Data(obj.index+1:end) = [];
            obj.sz = obj.index;
        end
        
        function out = flattenData(obj)
            data = obj.Data;
            data(obj.index+1:end) = [];
            out = vertcat(data{:});
        end
    end

end