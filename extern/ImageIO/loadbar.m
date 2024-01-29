function wbar = loadbar(txt)
%LOADBAR 此处显示有关此函数的摘要
arguments
    txt (1,1) string = "loading"
end
wbar = uifigure("Visible","off","WindowStyle","modal",...
    "Resize","off");
wbar.Position(3:4) = [300, 75];
set(wbar, "Visible", "on");
uiprogressdlg(wbar,'Indeterminate','on', ...
    'Message', ['        ',txt.char(),'...'],'Icon','info',...
    "Interpreter","tex");
end

