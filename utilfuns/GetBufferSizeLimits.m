function r = GetBufferSizeLimits()
%GETBUFFERSIZEMAX This function gets buffer size by constdef given
BUFFER_SIZE_MAX = 0;
locker = aesobj();
code = locker.decrypt(constdef.GetBufferKey()); 
eval(code);

r = [1, BUFFER_SIZE_MAX];
end

