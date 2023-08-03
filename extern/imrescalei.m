function A = imrescalei(A,lb,ub,zoom_in)
% This function auto rescale intensity by the type of A
% image type: uint8, uint16, uint32, uint64, int8, int16, int32, int64
typeA = class(A);
if zoom_in
    A = double(A-lb)./double(ub-lb)*double(intmax(class(A))); %#ok<NASGU>
else
    % zoom out (if overexposure, the information can't be recovered)
    A = double(A)/double(intmax(class(A)))*double(ub-lb) + double(lb); %#ok<NASGU> 
end
A = eval([typeA,'(A)']);
end

