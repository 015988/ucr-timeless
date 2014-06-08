clear
b=tic;



if (exist('intrinsic.mat','file'))
    disp('Loading intrinsic.mat')
    load('intrinsic')
else
    conf
    %flog=fopen('log.txt','w');
    flog=1;
    
    fprintf(flog,'starting intrinsic calibration\n');
    fprintf(flog,'%s\n',datestr(now));
    for cv=1:length(calVid)
        
        dat=load(calDat{cv});
        cI=zeros(2*size(dat,1),size(dat,2)/2);
        for i=1:size(dat,1)
            cI(2*(i-1)+1,:)=dat(i,1:2:end);
            cI(2*(i-1)+2,:)=dat(i,2:2:end);
        end
        ref=load(calRef{cv});
        cW=zeros(2*size(ref,1),size(ref,2)/3);
        for i=1:size(ref,1)
            cW(2*(i-1)+1,:)=ref(i,1:3:end);
            cW(2*(i-1)+2,:)=ref(i,2:3:end);
        end
        
        vCal=VideoReader(calVid{cv}); %#ok<TNMLP>
        [A(cv).A A(cv).kk]=IntrinsicCalibration(cI,cW,vCal.Height,vCal.Width,1);
        fprintf(flog,'Done %g\n',cv);
    end
    save -v7.3 intrinsic
end

%%
if (exist('extrinsic.mat','file'))
    disp('Loading extrinsic.mat\n')
    load('extrinsic.mat');
else
    conf
    %flog=fopen('log.txt','w');
    flog=1;
    fprintf(flog,'Extrinsic\n');
    fprintf(flog,'%s\n',datestr(now));
    
    for v=1:length(datName)
        
        dat=load(datName{v});
        cI=zeros(2*(vid(v).range(2)-vid(v).range(1)+1),size(dat,2)/2);
        ci=1;
        for i=vid(v).range(1):vid(v).range(2)
            cI(2*(ci-1)+1,:)=dat(i,1:2:end);
            cI(2*(ci-1)+2,:)=dat(i,2:2:end);
            ci=ci+1;
        end
        ref=load(refName{v});
        cW=zeros(2*(vid(v).range(2)-vid(v).range(1)+1),size(ref,2)/3);
        ci=1;
        for i=vid(v).range(1):vid(v).range(2)
            cW(2*(ci-1)+1,:)=ref(i,1:3:end);
            cW(2*(ci-1)+2,:)=ref(i,2:3:end);
            ci=ci+1;
        end
        
        clb(v).clb=PrepareCalibration(cI,cW,A(v).A,A(v).kk,vid(v).static);
        fprintf(flog,'Done %g\n',v);
    end
    save -v7.3 extrinsic
end

%%
conf
%flog=fopen('log.txt','w');
flog=1;
for v=1:length(datName)
    curStep=2*ceil(vidStepPerc*((1+vid(v).range(2)-vid(v).range(1))/100)/2);
    curRange(v).rangeMin=vid(v).range(1):round(curStep/2):vid(v).range(2);
    curRange(v).rangeMax=min(vid(v).range(2),curRange(v).rangeMin+curStep);
    data(v).v=VideoReader(vids{v}); %#ok<TNMLP>
end
%%

if 0 %debug purposes = to check orientation from both cameras
    for v=1:length(vids) %#ok<UNRCH>
        vr(1)=curRange(v).rangeMin(curBlock);
        vr(2)=curRange(v).rangeMax(curBlock);
        tmp=read(data(v).v,vr);
        for f=1:size(tmp,4)
            tPoints=round(project_points(clb(v).clb(vr(1)+f-curRange(v).rangeMin(1)).clb,[0 0 0;19.1 0 0;0 2*19.1 0]'));
            tPoints(tPoints<=0)=2;
            tPoints(1,tPoints(1,:)>size(tmp,2))=2;
            tPoints(2,tPoints(2,:)>size(tmp,1))=2;
            tImg=squeeze(tmp(:,:,:,f));
            tImg((tPoints(2,1)-1):(tPoints(2,1)+1),(tPoints(1,1)-1):(tPoints(1,1)+1),1)=255;
            tImg((tPoints(2,1)-1):(tPoints(2,1)+1),(tPoints(1,1)-1):(tPoints(1,1)+1),2)=0;
            tImg((tPoints(2,1)-1):(tPoints(2,1)+1),(tPoints(1,1)-1):(tPoints(1,1)+1),3)=0;
            
            tImg((tPoints(2,2)-1):(tPoints(2,2)+1),(tPoints(1,2)-1):(tPoints(1,2)+1),1)=0;
            tImg((tPoints(2,2)-1):(tPoints(2,2)+1),(tPoints(1,2)-1):(tPoints(1,2)+1),2)=255;
            tImg((tPoints(2,2)-1):(tPoints(2,2)+1),(tPoints(1,2)-1):(tPoints(1,2)+1),3)=0;
            
            tImg((tPoints(2,3)-1):(tPoints(2,3)+1),(tPoints(1,3)-1):(tPoints(1,3)+1),1)=0;
            tImg((tPoints(2,3)-1):(tPoints(2,3)+1),(tPoints(1,3)-1):(tPoints(1,3)+1),2)=0;
            tImg((tPoints(2,3)-1):(tPoints(2,3)+1),(tPoints(1,3)-1):(tPoints(1,3)+1),3)=255;
            imwrite(tImg,sprintf('ex%03d_%03d.png',v,vr(1)+f-1));
        end
    end
end
%%
cor(1)='k';
cor(2)='b';
allPoints=[];
for curBlock=1:length(curRange(refVid).rangeMin)
    fprintf(flog,'Loading videos - block %03d\n',curBlock);
    fprintf(flog,'%s\n',datestr(now));
    
    
    for v=1:length(vids)
        vr(1)=curRange(v).rangeMin(curBlock);
        vr(2)=curRange(v).rangeMax(curBlock);
        
        tmp=double(read(data(v).v,vr))/255.0;
        %tTarget=double(squeeze(tmp(:,:,1,:)));
        % BLUE -> tTarget=double(squeeze(tmp(:,:,3,:))-squeeze(tmp(:,:,1,:)));
        % RED -: tTarget=double(squeeze(tmp(:,:,1,:))-squeeze(tmp(:,:,2,:)));
        tTarget=double(squeeze(tmp(:,:,3,:))-squeeze(tmp(:,:,1,:)));
        tTarget(tTarget<0)=0;
        %tTarget=tTarget/max(max(max(tTarget)));
        tRed=double(squeeze(min(tmp,[],3)));
        tRed=1-tRed;%/max(max(max(tRed)));
        tTarget=tTarget.*tRed;
        tmp=permute(tmp,[1 2 4 3]);
        tmpSize=size(tmp);
        tmp=reshape(tmp,[numel(tmp)/3 3]);
        tmp=reshape(rgb2hsv(tmp),[tmpSize(1) tmpSize(2) tmpSize(3) 3]);
        tmp=squeeze(tmp(:,:,:,2));
        data(v).im=tmp.*tTarget;%/max(max(max(tmp))); %#ok<*SAGROW>
        data(v).im=data(v).im/max(max(max(data(v).im)));
        %figure,imshow(data(v).im),pause
        data(v).im=imdilate(data(v).im,ones(6,6,1));
        data(v).im=imclose(data(v).im,ones(6,6,3));
        data(v).im(data(v).im<(threshold/5))=0;
        %        for i=1:size(data(v).im,3)
        %            imwrite(squeeze(data(v).im(:,:,i)),sprintf('im%01d_%04d.png',v,vr(1)+i-1));
        %        end
        fprintf(flog,'Done %g',v);
        %figure,imshow(data(v).im(:,:,1),[]);
    end
    fprintf(flog,'%s\n',datestr(now));
    fprintf(flog,'Reconstructing points\n');
    
    for v=1:length(data)
        for f=1:size(data(v).im,3)
            %curClb(v).clb(f).clb=clb(v).clb(vr(1)+f-1).clb;
            curClb(v).clb(f).clb=clb(v).clb(curRange(v).rangeMin(curBlock)-curRange(v).rangeMin(1)+f).clb;
        end
    end
    
    %tmp=reconstructFromImages(curClbs,data,f,1000*[-2 0 -2 2 -1 1],0.75,150E3)';
    a=tic;
    recPoints=unsyncReconstructFromImages(curClb,data,1000*[-2 2 -2 2 -2 2],threshold,nPoints,minPoints,2,refVid,outRange)';
    toc(a)
    recPoints(:,5)=recPoints(:,5)+curRange(refVid).rangeMin(curBlock)-1;
    

    %figure(1),plot3(recPoints(:,1),recPoints(:,2),recPoints(:,5),'MarkerSize',1,'Marker','.','Color',cor(mod(curBlock,2)+1),'LineStyle','none')
    %hold on, grid on, axis tight square
    %figure(1),plot3(recPoints(recPoints(:,4)>0,1),recPoints(recPoints(:,4)>0,2),recPoints(recPoints(:,4)>0,5),'MarkerSize',1,'Marker','.','Color',cor(mod(curBlock,2)+1),'LineStyle','none')
    %hold on, grid on,
    %plot3(recPoints(recPoints(:,4)<=0,1),recPoints(recPoints(:,4)<=0,2),recPoints(recPoints(:,4)<=0,5),'MarkerSize',1,'Marker','.','Color','r','LineStyle','none')
    %axis tight square
    
    %save(sprintf('%03dto%03d.mat',curRange(refVid).rangeMin(curBlock),curRange(refVid).rangeMax(curBlock)),'recPoints','allPoints');    
    allPoints=[allPoints;recPoints]; %#ok<AGROW>
    
    %axis tight equal
    %clear data;
    fprintf(flog,'%s\n',datestr(now));
end
save bom
save('pontos.mat','allPoints');
datestr(now)
toc(b)
fprintf(flog,'%s\n',datestr(now));
fclose(flog);
%
% fprintf(flog,'labelling')
% kn=knnsearch(recPoints,recPoints,'K',5);
% lab=zeros(size(recPoints,1),1);
% curId=0;
% for i=1:length(lab)
%     if (lab(i)==0)
%         curId=curId+1;
%         lab(i)=curId;
%
%         Q=kn(i,2:end);
%         while (~isempty(Q))
%             c=Q(1);
%             Q=Q(2:end);
%             if (lab(c)==0)
%                 lab(c)=curId;
%                 Q=[Q kn(c,2:end)]; %#ok<AGROW>
%             end
%         end
%     end
% end
%
% %%
% fprintf(flog,'%s\n',datestr(now));
% fprintf(flog,'showing')
% figure(1), hold on,grid on
% figure(2), hold on,grid on
%
% for i=1:curId
%     if (sum(lab==i)>100)
%         figure(1),plot3(recPoints(lab==i,1),recPoints(lab==i,2),recPoints(lab==i,3),'+','MarkerEdgeColor',[rand(1,1) rand(1,1) rand(1,1)])
%         axis tight equal
%         cI=recPoints(lab==i,1); y = recPoints(lab==i,2); z= recPoints(lab==i,3);
%
%         k=convhull(cI,y,z);
%         k=unique(reshape(k,[1 numel(k)]));
%
%         cI=cI(k);y=y(k);z=z(k);
%
%         dt = DelaunayTri(cI,y,z);
%         Tes = dt(:,:);
%         cW = [cI(:) y(:) z(:)];
%         figure(2),tetramesh(Tes,cW,'FaceAlpha',1,'EdgeColor','none','FaceColor',[rand(1,1) rand(1,1) rand(1,1)]);
%         axis tight equal
%     else
%         lab(lab==i)=0;
%     end
% end
% axis tight equal, grid on
% %%
% figure
% wCam=2;
% for f=1:size(data(wCam).im,3)
%     imshow(data(wCam).im(:,:,f)),hold on
%     pointsChosen=recPoints(:,4)==f;
%     points=recPoints(pointsChosen,1:3);
%     nlab=lab(pointsChosen);
%     proj=round(project_points(curClbs(wCam).clb,points));
%     for i=1:curId
%         if (sum(lab==i)>100)
%             plot(proj(1,nlab==i),proj(2,nlab==i),'+','color',rand(1,3))
%         end
%     end
%     axis tight
%     pause
% end
%
%
