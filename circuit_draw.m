function [] = circuit_draw(cmap)

for i = 1:size(cmap,1)
    for j = i+1:size(cmap,1)
        assert(cmap(i,j)==cmap(j,i),'contact map is not symmetric')
    end
end

for i = 1:size(cmap,1)
    plot(i,0,'o','color','k','Linewidth',2)
    hold on
end
plot([1,size(cmap,1)],[0,0],'-k')

for i = 1:size(cmap,1)
    for j = (i + 1):size(cmap,1)
        if cmap(i,j) > 0
            X = [i mean([i,j]) j];
            Y = [0 0.1*(j - i) 0];
            t = linspace(1,3,30);
            xx = spline(1:3,X,t);
            yy = spline(1:3,Y,t);
            plot(xx,yy,'k','Linewidth',log(cmap(i,j)+1));
        end
    end
end

axis([1 size(cmap,1) -0.01*size(cmap,1) 0.1*size(cmap,1)])
set(gca,'Fontsize',20)
set(gca,'YTickLabel','')
xlabel('Node')
set(gca,'XTick',1:size(cmap,1));
end
