function u=generateProjection(A,R,T,W)
if (size(W,1)==3)
    W(4,:)=1;
end
u=A*[R T]*W;
u=u(1:2,:)./repmat(u(3,:),[2 1]);

%im=zeros(ceil(max(u,[],2)+20)');
%for i=1:size(u,2)
%    im(round(u(1,i)),round(u(2,i)))=1;
%end
