classdef weightClass < handle
    % takes weight and location of object and processes that information
    
    properties
        name
        weight
        location
        moment
    end
    
    methods
        function obj = weightClass(n,w,x,L)
            obj.name = n;
            if length(w) > 1
                obj.weight = sum(w);
            else
                obj.weight = w;
            end
            obj.location = x*L;
            obj.moment = w*x*L;
        end
        function [w, m] = totalWM(obj)
            w = sum([obj.weight]);
            m = sum([obj.moment]);
        end
        function lb2kg(obj)
            for i = 1:size(obj,2)
                obj(i).weight = obj(i).weight * 0.453592;
            end
        end
        function updateWeight(obj,weight)
            obj.weight = weight;
            obj.moment = weight * obj.location;
        end
    end
end
