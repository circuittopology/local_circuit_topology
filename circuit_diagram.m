% pdb file should be a single chain, with hydrogens removed
% stride file should correspond exactly to the pdb file
function [cmap3, s, segment, numbering] = circuit_diagram(pdb_file, stride_file, ss_elements, cutoff_distance, cutoff_numcontacts, plotcolor, plot_figs, plot_circle, plot_ions, plot_disulfide)
switch nargin
    case 2
        ss_elements = {'H','E','B','b','G'};
        %ss_elements = {'E','H'};
        cutoff_distance = 3.8; 
        cutoff_numcontacts = 6; 
        %plotcolor = [0.7 0.7 0.7];
        plotcolor = 'k';
        plot_figs = 1;
        plot_circle = 0;
        plot_ions = 0;
        plot_disulfide = 0;
end

        
% read pdb file
prot = pdbread(pdb_file);

% read Stride secondary structure file
fid = fopen(stride_file);
s = [];
seq = [];
while 1
    tline = fgets(fid);
    if strcmp(tline(1:3),'SEQ')
        seq = [seq, tline(11:60)];
    end
    if strcmp(tline(1:3),'STR')
        s = [s, tline(11:60)];
    end
    if strcmp(tline(1:3),'LOC')
        break
    end
end
seq = seq(isstrprop(seq,'alpha')); % like Python strip()
fclose(fid);
if length(s)>length(seq)
    s = s(1:length(seq));
end


% protein info
residue_number = [prot.Model.Atom.resSeq];
xcoor = [prot.Model.Atom.X];
ycoor = [prot.Model.Atom.Y];
zcoor = [prot.Model.Atom.Z];
name = {prot.Model.Atom.AtomName};
seq_length = length(s);


% number residues 1 to end, skipping over gaps.
numbering = unique(residue_number);
[residue_number,y] = find(residue_number==numbering'); 
assert(length(numbering)==seq_length, "different number of residues in stride file vs. pdb file")


% divide residues into segments based on secondary structure
nseg = 1;
segment = zeros(1,seq_length);
for i = 1:seq_length
    if ismember(s(i),ss_elements)
        segment(i) = nseg;
        if i == seq_length
            nseg = nseg + 1;
        elseif s(i + 1) ~= s(i)
            nseg = nseg + 1;
        end
    end
end
nseg = nseg - 1;

if plot_figs
    figure(1)
    imagesc(segment)
    colors = colormap(jet);
    colors(1,:) = [1,1,1];
    colormap(colors)
    set(gcf,'Position',[0 0 600 50])
    set(gca,'ytick',[])
    title('Sequence')
end


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
        if cmap(i,j) == 1
            seg_i = segment(residue_number(i));
            seg_j = segment(residue_number(j));
            if (seg_i ~= 0) && (seg_j ~= 0) 
                cmap2(seg_i,seg_j) = cmap2(seg_i,seg_j)+1;
            end
        end
    end
end

% make symmetric and set diagonal equal to zero
cmap2 = cmap2 + cmap2'; 
cmap2 = cmap2 - diag(diag(cmap2));

if plot_figs
    figure(2)
    imagesc(log(cmap2))
    axis square
    colormap bone
    title('log contact map before cutoff','Fontsize',15)
end

% circuit diagram based on segment-segment contact map
if plot_figs
figure(3)
hold on
% for i = 1:(nseg)
%    plot(mean(find(segment==i)),0,'o','color',plotcolor,'Linewidth',2,'MarkerFaceColor',[0.5 0.5 0.5])
%    hold on
% end
plot([1,seq_length],[0,0],'-k')

for i = 1:nseg
    for j = (i + 1):nseg
        if cmap2(i,j) >= cutoff_numcontacts
            X = [mean(find(segment==i)) mean(find(segment==j))];
            r = (X(2)-X(1))/2;
            center = mean(X);
            xx = X(1):0.01:X(2);
            yy = sqrt(r^2 - (xx - center).^2);
            plot(xx,yy,'color',plotcolor,'Linewidth',log(cmap2(i,j)+1));
        end
    end
end
%set(gca,'visible','off')
axis([0 length(s) 0 length(s)/2])
pbaspect([2 1 1])
set(gca,'Fontsize',20)
set(gca,'YTick','')
xlabel('Residue')
hold on
end

cmap3 = (cmap2 >= cutoff_numcontacts);



% ion information
if plot_ions
    ion_list = {'NA','K','FE','CA','CL','MG','ZN'};
    if isfield(prot.Model,"HeterogenAtom")
        ion_name = {prot.Model.HeterogenAtom.element};
        ion_xcoor = [prot.Model.HeterogenAtom.X];
        ion_ycoor = [prot.Model.HeterogenAtom.Y];
        ion_zcoor = [prot.Model.HeterogenAtom.Z];
        
        % ions only
        ision = ismember(ion_name, ion_list);
        ion_name = ion_name(ision==1);
        ion_xcoor = ion_xcoor(ision==1);
        ion_ycoor = ion_ycoor(ision==1);
        ion_zcoor = ion_zcoor(ision==1);
        
        % we'll allow contact with any partial occupancy location for each ion
        ion_number = [prot.Model.HeterogenAtom.resSeq];
        ion_numbering = unique(ion_number);
        [ion_number,y] = find(ion_number==ion_numbering');
        
        ion_cmap = pdist2([xcoor; ycoor; zcoor]',[ion_xcoor;ion_ycoor;ion_zcoor]') < cutoff_distance;
        
        % residue-ion contact map, based on atom-atom contact map
        ion_cmap2 = zeros(seq_length,max(ion_number));
        for i = 1:natoms
            for j = 1:length(ion_name)
                if ion_cmap(i,j) == 1
                    seg_i = residue_number(i);
                    seg_j = ion_number(j);
                    if (seg_i ~= 0) && (seg_j ~= 0)
                        ion_cmap2(seg_i,seg_j) = ion_cmap2(seg_i,seg_j)+1;
                    end
                end
            end
        end
        
        ion_cutoff_numcontacts = 1;
        for nion = 1:max(ion_number)
            for i = 1:seq_length
                for j = (i + 1):seq_length
                    if (ion_cmap2(i,nion) >= ion_cutoff_numcontacts) && (ion_cmap2(j,nion) >= ion_cutoff_numcontacts)
                        X = [i j];
                        r = (X(2)-X(1))/2;
                        center = mean(X);
                        xx = X(1):0.01:X(2);
                        yy = sqrt(r^2 - (xx - center).^2);
                        plot(xx,-yy,'g','Linewidth',log((ion_cmap2(i,nion) + ion_cmap2(j,nion))+1));
                    end
                end
            end
        end
        axis([0 length(s) -length(s)/2 length(s)/2])
        pbaspect([1 1 1])
        
    end
end


% disulfide bonds
if plot_disulfide
    ds_cutoff = 2.8;
    
    allx = [prot.Model.Atom.X];
    ally = [prot.Model.Atom.Y];
    allz = [prot.Model.Atom.Z];
    
    iss = strcmp({prot.Model.Atom.AtomName},'SG').*strcmp({prot.Model.Atom.resName},'CYS');
    scoor = [allx(iss==1); ally(iss==1); allz(iss==1)];
    s_resnum = residue_number(iss==1);
    duplicate_s = zeros(1,length(s_resnum));
    for i = 2:length(duplicate_s)
        if s_resnum(i)==s_resnum(i-1)
            duplicate_s(i) = 1;
        end
    end
    s_resnum(duplicate_s==1) = [];
    scoor(:,duplicate_s==1) = [];
    
    ds_cmap = squareform(pdist(scoor')) < ds_cutoff;
    
    seq_length
    ds_cmap2 = zeros(seq_length);
    for i = 1:length(s_resnum)
        for j = 1:length(s_resnum)
            if ds_cmap(i,j) == 1
                seg_i = s_resnum(i);
                seg_j = s_resnum(j);
                if (seg_i ~= 0) && (seg_j ~= 0)
                    [seg_i, seg_j]
                    ds_cmap2(seg_i,seg_j) = ds_cmap2(seg_i,seg_j)+1;
                end
            end
        end
    end
    
        for i = 1:seq_length
            for j = (i + 1):seq_length
                if ds_cmap2(i,j) > 0
                    X = [i j];
                    r = (X(2)-X(1))/2;
                    center = mean(X);
                    xx = X(1):0.01:X(2);
                    yy = sqrt(r^2 - (xx - center).^2);
                    plot(xx,-yy,'color',[0.5 0.5 0.5],'Linewidth',log(4));
                end
            end
        end
            
    axis([0 length(s) -length(s)/2 length(s)/2])
    pbaspect([1 1 1])
    
end



% circular diagram
if plot_circle
    figure(4)
    plot(cos(0:pi/50:2*pi),sin(0:pi/50:2*pi),'color',[0.8 0.8 0.8])
    hold on
    for i = 1:nseg
        for j = (i + 1):nseg
            if cmap2(i,j) >= cutoff_numcontacts
                theta1 = mean(find(segment==i))*2*pi()/seq_length;
                theta2 = mean(find(segment==j))*2*pi()/seq_length;
                X = [cos(theta1),cos(theta2)];
                Y = [sin(theta1),sin(theta2)];
                plot(X,Y,plotcolor,'Linewidth',log(cmap2(i,j)+1));
                hold on
            end
        end
    end
    symbol = containers.Map;
    symbol('H') = 'o';
    symbol('E') = '^';
    symbol('B') = '*';
    symbol('b') = '*';
    symbol('G') = '*';
    keys = symbol.keys;
    for i = 1:length(ss_elements)
        assert(ismember(ss_elements{i},keys), "Add a symbol for missing element in the map for the circle diagram")
    end
    colors2 = jet(seq_length);
    for i = 1:nseg
        theta = mean(find(segment==i))*2*pi()/seq_length;
        plot(cos(theta),sin(theta),symbol(s(find(segment==i,1))),'color',colors2(floor(mean(find(segment==i))),:),'Linewidth',2,'Markersize',15)
        hold on
    end
    axis square
    axis([-1 1 -1 1])
    set(gca,'xtick','')
    set(gca,'ytick','')
end