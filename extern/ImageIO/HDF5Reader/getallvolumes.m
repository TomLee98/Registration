function [img,time]=getallvolumes(fname,lasers,piezo,img_idx,res,volumes,t,img)
    
        
        [~,~,int_piezo]=unique(piezo);
        %volumes=volumes(img_idx);
        %piezo = piezo(img_idx);
        %lasers = lasers(img_idx);
        use_frames=volumes~=0;
        %frames=find(volumes==volume);
       
        %offset=prod(res)*(frames-1)+1;
        %start=[1,offset(1)];
        %count=[1,prod(res)*length(];
        
        data=h5read(fname,'/img');

        unique_lasers=find(any(lasers,1));
        
        data3d=reshape(data,res(1),res(2),length(img_idx));
        
        which_lasers=zeros(size(lasers,1),1);
        for ii=1:length(lasers)
            tmp=find(lasers(ii,:));
            if ~isempty(tmp)
                which_lasers(ii)=find(lasers(ii,:));
            end
        end
        
        [~,~,which_laser_idx]=unique(which_lasers(use_frames));
        data3d_use=data3d(:,:,use_frames);
        int_piezo_use=int_piezo(use_frames);
        volumes=volumes(use_frames);
        for ii=1:length(which_laser_idx)
%             img(:,:,int_piezo_use(ii),volumes(ii),which_laser_idx(ii))=data3d_use(:,:,ii); %xyztc
            img(:,:,int_piezo_use(ii),which_laser_idx(ii),volumes(ii))=data3d_use(:,:,ii); %xyzct
        end
        unique_volumes=unique(volumes(volumes>0));
        for ii=1:length(unique_volumes)
           time(ii)=mean(t(volumes==unique_volumes(ii)));
            
        end
        img = flip(flip(permute(img,[2,1,3,4,5]),1),3);
        
    