function [mat, c] = get_matrix(cmap,plot_figs)
% contacts are numbered left to right. 
% if two contacts begin at the same point, the one that ends first comes
% first.
switch nargin
    case 1
        plot_figs = 1;
end

assert(trace(cmap) == 0, "there should be zeros on the diagonal")
[rows, cols] = find(triu(cmap));
[rows,idx] = sort(rows);
cols = cols(idx);
mat = eye(length(rows))*8;
c = cell(length(rows));
for z = 1:length(rows)
    c{z,z} = 'P';
end
for x = 1:length(rows)
    for y = (x+1):length(rows)
        k = rows(x);
        m = cols(x);
        i = rows(y);
        j = cols(y);
        if k < i && m > j
            mat(x,y) = 2;
            mat(y,x) = 1;
            c{x,y} = 'P-1';
            c{y,x} = 'P';
        elseif m < i
            mat(x,y) = 3;
            mat(y,x) = 3;
            c{x,y} = 'S';
            c{y,x} = 'S';
        elseif (k == i && m < j) 
            mat(x,y) = 5;
            mat(y,x) = 6;
            c{x,y} = 'CP';
            c{y,x} = 'CP-1';
        elseif k == i && m > j || (m == j)
            mat(x,y) = 6;
            mat(y,x) = 5;
            c{x,y} = 'CP-1';
            c{y,x} = 'CP';
        elseif m == i
            mat(x,y) = 7;
            mat(y,x) = 7;
            c{x,y} = 'CS';
            c{y,x} = 'CS';
        else
            assert(m > i && j > m,'something went wrong')
            mat(x,y) = 4;
            mat(y,x) = 4;
            c{x,y} = 'X';
            c{y,x} = 'X';
        end
    end
end

if plot_figs
    figure
    circuit_draw(cmap)
    for x = 1:length(rows)
        text((rows(x)+cols(x))/2,0.1*(cols(x)-rows(x)),sprintf('%g',x),'Fontsize',25)
    end
    figure
    imagesc(mat)
    caxis([1 8])
    colormap([172/255 200/255 247/255; 174/255 213/255 129/255; 131/255 139/255 197/255; 186/255 155/255 201/255; 172/255 200/255 247/255; 174/255 213/255 129/255; 131/255 139/255 197/255; 218/255 219/255 228/255]);
    
    % P, P-1, S, X, CP, CP-1, CS
    axis square
    set(gca,'xtick',1:1:length(rows))
    set(gca,'ytick',1:1:length(rows))
    set(gca,'Fontsize',20)
    
    hold on
    y = 0.5;
    ymax = length(rows)+1;
    while y < ymax
        %plot([0 ymax],[y y],'k','Linewidth',1.5)
        %plot([y y],[0 ymax],'k','Linewidth',1.5)
        y = y + 1;
    end
    
%     figure
%     imagesc([3 1 2 4]')
%     caxis([1 8])
%     axis square
%     pbaspect([1 4 1])
%     set(gca,'visible','off')
%     colormap([172/255 200/255 247/255; 174/255 213/255 129/255; 131/255 139/255 197/255; 186/255 155/255 201/255; 174/255 213/255 129/255; 174/255 213/255 129/255; 131/255 139/255 197/255; 218/255 219/255 228/255]);
%     hold on
%     y = 0.5;
%     while y < 5
%         plot([0 5],[y y],'k','Linewidth',1.5)
%         y = y + 1;
%     end
%     plot([0.5 0.5],[0 5],'k','Linewidth',1.5)
%     plot([1.5 1.5],[0 5],'k','Linewidth',1.5)
end
    figure(2)