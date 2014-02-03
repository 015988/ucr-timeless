clear
delete log.txt

nPoints=2000;
epRange=[10 5 2 1 0.5];

pStep=20;
minPoints=500;

minV=[ 1  1  1]'-1;
maxV=[11 11 11]'+1;


rng(9)
W=generate3d(6,300,minV'+1,maxV'-1);
figure(9),plot3(W(1,:),W(2,:),W(3,:),'.'),hold on, grid on, axis tight equal
figure(10),plot(sqrt(sum((W(:,1:(size(W,2)-1))-W(:,2:end)).^2))*100,'+-'),hold on, grid on %speed

A1=[100 0 400;0 100 800;0 0 1];
% T1=[-1 8 1/2]'; standard
% T1=[1 1 8]'; small baseline
% T1=[8 1 4]'; 
T1=-[30 5 20]'; 
R1=eye(3); 
u1=generateProjection(A1,R1,T1,W(:,1:3:end));
w1size=((size(u1,2)*(pStep)/100)/2);
w1start=1:floor(w1size):size(u1,2);
w1end=w1start+2*ceil(w1size);
w1end(w1end>size(u1,2))=size(u1,2);
figure,plot(u1(1,:),u1(2,:),'+-');axis equal


A2=[100 0 400;0 100 800;0 0 1];
R2=eye(3);  %T2=[0 0 8]'; % T2=[8 -1 1]';
%T2=[8 1 -4]'; 
T2=-[5 30 20]'; 
u2=generateProjection(A2,R2,T2,W(:,2:3:end));
w2size=(size(u2,2)*(pStep)/100)/2;
w2start=1:floor(w2size):size(u2,2);
w2end=w2start+2*ceil(w2size);
w2end(w2end>size(u2,2))=size(u2,2);
figure,plot(u2(1,:),u2(2,:),'+-');

A3=[200 0 200;0 200 200;0 0 1];
R3=eye(3); 
T3=-[-20 -20 20]';
u3=generateProjection(A3,R3,T3,W(:,3:3:end));
w3size=(size(u3,2)*(pStep)/100)/2;
w3start=1:floor(w3size):size(u3,2);
w3end=w3start+2*ceil(w3size);
w3end(w3end>size(u3,2))=size(u3,2);
figure,plot(u3(1,:),u3(2,:),'+-');


pRes=[];
tRes=[];

for s=1:size(w1start,2)
    mklog(sprintf('Starting block %g/%g',s,size(w1start,2)));
    
    u1n=u1(:,w1start(s):w1end(s));
    u2n=u2(:,w2start(s):w2end(s));
    u3n=u3(:,w3start(s):w3end(s));
    
    
    c1=csapi(1:size(u1n,2),u1n);
    c2=csapi(1:size(u2n,2),u2n);
    c3=csapi(1:size(u3n,2),u3n);
    
    points=[];
    
    for epsilon=epRange
        mklog(sprintf('--Starting epsilon = %g',epsilon));
        ok=0;
        repeats=0;
        pInRegion=-1;

        while (ok==0)
            ok=1;
            if (isempty(points))
                mklog(sprintf('Creating %u new random points (0) - repeat %g',(100+repeats)*nPoints,repeats));
                points=newPoints((10+repeats)*nPoints,minV,maxV);
                justNew=1;
            end
            
              
            while ((pInRegion~=0)&&(pInRegion<minPoints))
                if (justNew~=1)
                    if (pInRegion>0) %already made a pass and didn't find enough points
                        points=points(:,used);
                    end
                    if (size(points,2)>(nPoints/4))
                        mklog(sprintf('Reducing to %g (%g)',(nPoints/4),size(points,2)));
                        points=points(:,randperm(size(points,2),nPoints/4));
                    end
                    
                    if (size(points,2)<nPoints)
                        mklog(sprintf('Creating points in lines (%g)',size(points,2)));
                        points=[points pointsInLines(points)];
                        mklog(sprintf('Creating %g random points (%g)',2*nPoints-size(points,2),size(points,2)));
                        mP=mean(points,2);
                        sP=std(points,0,2);
                        sP(sP<0.1)=0.1;
                        points=[points normrnd(points,0.05) newPoints(nPoints-size(points,2),max(mP-10*sP,minV),min(mP+10*sP,maxV))];
                    end
                end
                justNew=0;
                mklog(sprintf('projections (%g)',size(points,2)));
                u1i=project(A1,R1,T1,points);
                mE=3+repeats;
                out1=(u1i(1,:)>max(u1n(1,:)+mE*epsilon)) | (u1i(2,:)>max(u1n(2,:)+mE*epsilon))|...
                    (u1i(1,:)<min(u1n(1,:)-mE*epsilon)) | (u1i(2,:)<min(u1n(2,:)-mE*epsilon));
                
                u2i=project(A2,R2,T2,points);
                out2=(u2i(1,:)>max(u2n(1,:)+mE*epsilon)) | (u2i(2,:)>max(u2n(2,:)+mE*epsilon))|...
                    (u2i(1,:)<min(u2n(1,:)-mE*epsilon)) | (u2i(2,:)<min(u2n(2,:)-mE*epsilon));
                
                  u3i=project(A3,R3,T3,points);
                  out3=(u3i(1,:)>max(u3n(1,:)+mE*epsilon)) | (u3i(2,:)>max(u3n(2,:)+mE*epsilon))|...
                      (u3i(1,:)<min(u3n(1,:)-mE*epsilon)) | (u3i(2,:)<min(u3n(2,:)-mE*epsilon));
                outAll= out1 | out2 | out3;
%                outAll= out1 | out2 ;
                used=find(outAll==0);
                pInRegion=length(used);
                
                figure(2),hold off
                plot3(points(1,outAll==0),points(2,outAll==0),points(3,outAll==0),'.b');
                hold on, grid on, axis tight equal
                plot3(W(1,:),W(2,:),W(3,:),'xk-','LineWidth',5);
                plot3(points(1,:),points(2,:),points(3,:),'or');
                mklog(sprintf('We have %g points in the region',pInRegion));
                print('-dpng','-r600',sprintf('zoom_%03.0f_%05.2f.png',s,epsilon));
                
                figure(3),hold off
                plot3(W(1,:),W(2,:),W(3,:),'.-k');
                hold on, grid on, axis tight equal
                plot3(points(1,out1==0),points(2,out1==0),points(3,out1==0),'.r');
                plot3(points(1,out2==0),points(2,out2==0),points(3,out2==0),'.g');
%                 plot3(points(1,out3==0),points(2,out3==0),points(3,out3==0),'.b');
                %print('-dpng','-r300',sprintf('cross_%03.0f_%05.2f.png',s,epsilon));
            end
            pInRegion=-1;
            mklog(sprintf('distance calculations (%g)',sum(outAll==0)));
            [nU1 near1]=pointsToCurveDist(u1i,c1,epsilon,find(outAll==0),epsilon==epRange(end));
            nU2=pointsToCurveDist(u2i,c2,epsilon,nU1,epsilon==epRange(end));
            nU3=pointsToCurveDist(u3i,c3,epsilon,intersect(nU1,nU2),epsilon==epRange(end));
            
            %used=nU2;
             used=nU3;
            figure(1),hold off,
            plot3(W(1,:),W(2,:),W(3,:),'-xk'),hold on, grid on
            plot3(points(1,used),points(2,used),points(3,used),'ro');
            axis tight equal
            %print('-dpng','-r300',sprintf('partial_%03.0f_%05.2f.png',s,epsilon));
%            plot3(points(1,:),points(2,:),points(3,:),'.b');
%            print('-dpng','-r300',sprintf('debug_%03.0f_%05.2f.png',s,epsilon));
%                print('-deps',sprintf('partial_%03.0f_%05.2f.eps',s,epsilon));
            print('-dpng','-r600',sprintf('partial_%03.0f_%05.2f.png',s,epsilon));
            
            if (epsilon==epRange(end)) %we only need the time at the end..
                mklog(sprintf('Computing timeRef (%g)',length(used)));
                timeRef=zeros(size(used,2),1);
                ii=1;
                for p=used
                    timeRef(ii)=near1{find(nU1==p,1)}+(w1start(s))-1;
                    ii=ii+1;
                end
            end
            if (length(used)<minPoints)
                ok=0;
                repeats=repeats+1;
            end
            points=points(:,used);
        end
    end
            
    mklog('block done');
    pRes=[pRes points];
    tRes=[tRes; timeRef];
    
    figure(4),hold off
    plot3(W(1,:),W(2,:),W(3,:),'xk-');
    hold on, grid on, axis tight equal
    plot3(pRes(1,:),pRes(2,:),pRes(3,:),'.');
    %print('-dpng','-r300',sprintf('total_%3.0f.png',s));
end
       
save preNC
spRes=csaps(tRes,pRes,0.025);
uRes=fnval(spRes,1:size(u1,2));
figure,plot3(uRes(1,:),uRes(2,:),uRes(3,:),'r+-')
hold on, grid on
plot3(W(1,:),W(2,:),W(3,:),'-+k'),hold on, grid on
axis tight equal
figure,plot(sqrt(sum((W(:,1:3:end)-uRes).^2))),hold on, grid on
save posNC            
