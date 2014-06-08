function chooseOriginRef(vidName,datName,refName)
dat=load(datName);
nFrames=size(dat,1)/2;
x=zeros(nFrames,size(dat,2)/2,2);
for i=1:size(dat,1)
    x(i,:,1)=dat(i,1:2:end);
    x(i,:,2)=dat(i,2:2:end);
end

% for i=1:(nFrames-1)
%     if (~(any(any(x(i,:,:)==-1)))&&...
%         ~(any(any(x(i+1,:,:)==-1))))
%         D=pdist2(squeeze(x(i,:,:)),squeeze(x(i+1,:,:)));
%         if (~all(all(eye(size(D,1))==(D==repmat(min(D,[],2),[1 size(D,1)])))))
%             dir
%         end
%     end
% end



ref=load(refName);
X=zeros(nFrames,size(ref,2)/3,2);
for i=1:size(ref,1)
    X(i,:,1)=ref(i,1:3:end);
    X(i,:,2)=ref(i,2:3:end);
end


v=VideoReader(vidName);

while(any(any(x(i,:,:)==-1))); 
    i=i+1;
end

figure,imshow(read(v,i));
hold on
for p=1:size(x,2)
    text(x(i,p,1),x(i,p,2),sprintf('%d',p),'Color','red');
end
pO=find(all((squeeze(X(i,:,:))==0)'));
plot(x(i,pO,1),x(i,pO,2),'ob');

pMaxX=find((X(i,:,1)==0)&(max(X(i,:,2))==X(i,:,2)));
pMaxY=find((X(i,:,2)==0)&(max(X(i,:,1))==X(i,:,1)));
plot(x(i,pMaxX,1),x(i,pMaxX,2),'xb');
plot(x(i,pMaxY,1),x(i,pMaxY,2),'+b');


end