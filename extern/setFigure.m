function setFigure(f,arg,scr)
%setFigure - 设置figure在屏幕中的位置和大小
%
%setFigure(f,arg)
%f: figure handle
%arg: [centerX,centerY,width,height], 使用归一化参数，center相对于屏幕中央为原点
%
%example: 以屏幕中心为原点，宽度占比100%，高度占比100%
%f = figure(1);
%setFigure(f,[0,0,1,1]) 

arguments
    f (1,1);
    arg (1,4) double {mustBeInRange(arg, 0, 2)};
    scr (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
end


screenSize = get(0,'ScreenSize');

%最大化窗口
if arg(3)==arg(4) && arg(3)==1
    warning('off');	            % 关闭相关的警告提示（因为调用了非公开接口）
    if scr ~= 0
        % move the figure to the place
        pos_x = screenSize(3)-screenSize(1)/2;
        pos_y = (screenSize(4)-screenSize(1))/2; % screen layout is [X,Y]
        set(f,'position',[pos_x,pos_y,1,1]);
    end
    jFrame = get(f,'JavaFrame');%#ok<JAVFM> % 获取底层 Java 结构相关句柄
    pause(0.5);					% 在 Win 10，Matlab 2017b 环境下不加停顿会报 Java 底层错误。各人根据需要可以进行实验验证
    set(jFrame,'Maximized',1);	%设置其最大化为真（0 为假）
    pause(0.5);					% 个人实践中发现如果不停顿，窗口可能来不及变化，所获取的窗口大小还是原来的尺寸。各人根据需要可以进行实验验证
    warning('on');		        % 打开相关警告设置
else
    width = arg(3)*(screenSize(3)-screenSize(1)+1);
    height = arg(4)*(screenSize(4)-screenSize(2)+1);
    pos_x = (screenSize(3)-screenSize(1)-width)/2 + (screenSize(3)-screenSize(1)+1)*scr;
    pos_y = (screenSize(4)-screenSize(1)-height)/2; % screen layout is [X,Y]
    set(f,'position',[pos_x,pos_y,width,height]);
end
end

