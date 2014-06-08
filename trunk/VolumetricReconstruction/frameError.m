function ferr=frameError(cI, cW, curIntrinsic, curOT)
if (size(curOT,1)>1)
    curOT=curOT';
end
%RT=[rodrigues(curOT(1:3)) curOT(4:6)'];
%curA=[curIntrinsic(1) 0 curIntrinsic(3);...
%    0 curIntrinsic(2) curIntrinsic(4);...
%    0 0 1];
tclb.KK=[curIntrinsic(1) 0 curIntrinsic(3);...
    0 curIntrinsic(2) curIntrinsic(4);...
    0 0 1];
tclb.RT=[rodrigues(curOT(1:3)) curOT(4:6)'];
tclb.dir.k=[curIntrinsic(5) curIntrinsic(6)];
tclb.dir.p=[0 0];
tclb.dir.s=[0 0];
xl=project_points(tclb,cW);

%         X=RT*[cW;zeros(1,size(cW,2));ones(1,size(cW,2))];
%         X=X(1:2,:)./repmat(X(3,:),[2 1]);
%         xl=X;
%
%         if ((size(xl,2)~=size(X,2))||(size(xl,1)~=2))
%             xl=NaN;
%         else
%             xl=curA*[xl;ones(1,size(xl,2))];
%             xl=xl(1:2,:)./repmat(xl(3,:),[2 1]);
%         end
%ferr=sum(sqrt(sum((xl-cI).^2)));
ferr=sum((sum((xl-cI).^2)));
global DEBUG_PLOT;
if (exist('DEBUG_PLOT','var'))
    if (DEBUG_PLOT==1)
        %figure(10),hold off,
        figure
        plot(cI(1,:),cI(2,:),'o');
        hold on, grid on
        plot(xl(1,:),xl(2,:),'r+');
        legend('Original value','Reprojected');
        title(sprintf('Error : %g',ferr));
    end
end

end
