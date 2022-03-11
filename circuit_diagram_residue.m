% PDB file should be a single chain, with hydrogens removed
function [cmap3, numbering] = circuit_diagram_residue(pdb_file, cutoff_distance, cutoff_numcontacts, plot_figs)
switch nargin
    case 1
        cutoff_distance = 5;
        cutoff_numcontacts = 6;
        plot_figs = 1;
end

prot = pdbread(pdb_file);


residue_number = [prot.Model.Atom.resSeq];
xcoor = [prot.Model.Atom.X];
ycoor = [prot.Model.Atom.Y];
zcoor = [prot.Model.Atom.Z];
name = {prot.Model.Atom.AtomName};
numbering = unique(residue_number);
[residue_number,y] = find(residue_number==numbering');



% divide into segments based on residue
nseg = length(numbering);
segment = 1:nseg;

% delete duplicate atoms due to multiple occupancy
duplicate = zeros(length(name));
for i = 2:length(duplicate)
    if (residue_number(i)==residue_number(i-1)) && strcmp(name(i),name(i-1))
        duplicate(i) = 1;
    end
end

residue_number(duplicate==1) = [];
xcoor(duplicate==1) = [];
ycoor(duplicate==1) = [];
zcoor(duplicate==1) = [];

natoms = length(residue_number);


% atom-atom contact map
cmap = squareform(pdist([xcoor; ycoor; zcoor]')) < cutoff_distance;

% segment-segment contact map, based on atom-atom contact map
cmap2 = zeros(nseg);
for i = 1:natoms
    for j = (i+1):natoms
        seg_i = segment(residue_number(i));
        seg_j = segment(residue_number(j));
        if cmap(i,j) == 1
            if (seg_i ~= 0) && (seg_j ~= 0)
                cmap2(seg_i,seg_j) = cmap2(seg_i,seg_j)+1;
            end
        end
    end
end

cmap2 = cmap2 + cmap2';

% set values close to diagonal to zero
for i = 1:length(cmap2)
    cmap2(i,i) = 0;
end
for i = 1:(length(cmap2)-1)
    cmap2(i,i+1) = 0;
    cmap2(i+1,i) = 0;
end
for i = 1:(length(cmap2)-2)
    cmap2(i,i+2) = 0;
    cmap2(i+2,i) = 0;
end
for i = 1:(length(cmap2)-3)
    cmap2(i,i+3) = 0;
    cmap2(i+3,i) = 0;
end

if plot_figs
    figure(2)
    imagesc(log(cmap2)) % set diagonal equal to zero
    axis square
    colormap bone
    title('log contact map before cutoff','Fontsize',15)
    
    
    % circuit diagram based on segment-segment contact map
    figure(3)
    %for i = 1:nseg
    %    plot(mean(find(segment==i)),0,'o','color',plotcolor,'Linewidth',2)
    %    hold on
    %end
    plot([1,length(numbering)],[0,0],'-k')
    hold on
    
    for i = 1:nseg
        for j = (i + 1):nseg
            if cmap2(i,j) >= cutoff_numcontacts
                X = [mean(find(segment==i)) mean(find(segment==j))];
                r = (X(2)-X(1))/2;
                center = mean(X);
                xx = X(1):0.01:X(2);
                yy = sqrt(r^2 - (xx - center).^2);
                plot(xx,yy,'k','Linewidth',log(cmap2(i,j)+1));
            end
        end
    end
    
    b = -0.25;
    set(gca,'Fontsize',20)
    set(gca,'YTickLabel','')
    xlabel('Residue')
    hold on
end

cmap3 = cmap2>=cutoff_numcontacts;



