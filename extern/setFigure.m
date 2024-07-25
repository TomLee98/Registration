function setFigure(f,arg,scr)
%setFigure - ����figure����Ļ�е�λ�úʹ�С
%
%setFigure(f,arg)
%f: figure handle
%arg: [centerX,centerY,width,height], ʹ�ù�һ��������center�������Ļ����Ϊԭ��
%
%example: ����Ļ����Ϊԭ�㣬���ռ��100%���߶�ռ��100%
%f = figure(1);
%setFigure(f,[0,0,1,1]) 

arguments
    f (1,1);
    arg (1,4) double {mustBeInRange(arg, 0, 2)};
    scr (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
end


screenSize = get(0,'ScreenSize');

%��󻯴���
if arg(3)==arg(4) && arg(3)==1
    warning('off');	            % �ر���صľ�����ʾ����Ϊ�����˷ǹ����ӿڣ�
    if scr ~= 0
        % move the figure to the place
        pos_x = screenSize(3)-screenSize(1)/2;
        pos_y = (screenSize(4)-screenSize(1))/2; % screen layout is [X,Y]
        set(f,'position',[pos_x,pos_y,1,1]);
    end
    jFrame = get(f,'JavaFrame');%#ok<JAVFM> % ��ȡ�ײ� Java �ṹ��ؾ��
    pause(0.5);					% �� Win 10��Matlab 2017b �����²���ͣ�ٻᱨ Java �ײ���󡣸��˸�����Ҫ���Խ���ʵ����֤
    set(jFrame,'Maximized',1);	%���������Ϊ�棨0 Ϊ�٣�
    pause(0.5);					% ����ʵ���з��������ͣ�٣����ڿ����������仯������ȡ�Ĵ��ڴ�С����ԭ���ĳߴ硣���˸�����Ҫ���Խ���ʵ����֤
    warning('on');		        % ����ؾ�������
else
    width = arg(3)*(screenSize(3)-screenSize(1)+1);
    height = arg(4)*(screenSize(4)-screenSize(2)+1);
    pos_x = (screenSize(3)-screenSize(1)-width)/2 + (screenSize(3)-screenSize(1)+1)*scr;
    pos_y = (screenSize(4)-screenSize(1)-height)/2; % screen layout is [X,Y]
    set(f,'position',[pos_x,pos_y,width,height]);
end
end

