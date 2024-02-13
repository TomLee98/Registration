classdef ImageWriter < handle
    %IMAGEWRITER This class is an image writer defination, which  support
    %mutiple image format writting
    
    properties
        Property1
    end
    
    methods
        function obj = ImageWriter(inputArg1,inputArg2)
            %IMAGEWRITER 构造此类的实例
            %   此处显示详细说明
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            outputArg = obj.Property1 + inputArg;
        end
    end
end

