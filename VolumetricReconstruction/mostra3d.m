function mostra3d(dados,cor,alfa)
%figure
try
    cor;
catch
    cor='red';
end

try 
    alfa;
catch
    alfa=0.2;
end
hold on
p1 = patch(isosurface(dados,.4),'FaceColor',cor,'EdgeColor','none','FaceAlpha',alfa);
p2= patch(isocaps(dados,.4),'FaceColor',cor,'EdgeColor','none','AlphaDataMapping','none','FaceAlpha',alfa);
isonormals(dados,p1);
view(3); axis vis3d tight
camlight('left','infinite');
camlight('right','infinite');
lighting gouraud
grid on
tam=size(dados);
axis([0 tam(2) 0 tam(1) 0 tam(3)]);
daspect([1 1 1])