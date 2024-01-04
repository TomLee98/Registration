function [img]=initializeimgs(volumes,lasers,res,img_idx)
% initializeimgs, determine the size of the img stack

    % x,y,z,c,t
    unique_lasers=find(any(lasers,1));
    num_lasers=length(find(any(lasers,1)));    
    
    for ii=1:size(lasers,2) 
        zdepth(ii)=length(find(lasers(volumes==1,ii)));
    end
    
%     img=int16(zeros(res(1),res(2),max(zdepth),length(unique(volumes(volumes>0))),num_lasers)); %xyztc
    img=int16(zeros(res(1),res(2),max(zdepth),num_lasers,length(unique(volumes(volumes>0))))); %xyzct
   