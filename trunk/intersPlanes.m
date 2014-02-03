function intersPlanes(W,C,edgeColor)
if (nargin==2)
    edgeColor='k';
end
    hold on,grid on
    D=W-repmat(C,[1 size(W,2)]);
    for i=2:size(W,2)

        patch([C(1) W(1,i-1)+D(1,i-1) W(1,i)+D(1,i) C(1)],...
              [C(2) W(2,i-1)+D(2,i-1) W(2,i)+D(2,i) C(2)],...
              [C(3) W(3,i-1)+D(3,i-1) W(3,i)+D(3,i) C(3)],...
              rand(4,1),...
              'FaceColor','none','EdgeColor',edgeColor);
    end
    plot3(W(1,:),W(2,:),W(3,:),'o-k','LineWidth',3);
end