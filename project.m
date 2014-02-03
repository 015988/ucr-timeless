function u=project(A,R,T,X)
P=A*[R T];
X(4,:)=1;
u=P*X;
u=u(1:2,:)./repmat(u(3,:),[2 1]);