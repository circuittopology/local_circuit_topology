files = dir('./out/cmap*');

for i = 1:length(files)
    i
    cmap3 = load(sprintf('./out/%s',files(i).name));
    c = get_matrix(cmap3,0);
    [rows, cols] = find(triu(cmap3));
    [rows,idx] = sort(rows);
    cols = cols(idx);
    
    % protein length and contact order
    %co = mean(abs(rows-cols));
    nres = length(cmap3);
    coi = zeros(1,length(cmap3));
    for j = 1:length(coi)
        if  ~isempty(mean(abs(find(cmap3(j,:))-j)))
            coi(j) = mean(abs(find(cmap3(j,:))-j));
        end
    end
    
    % local circuit topology
    p_number = zeros(1,length(cmap3));
    ip_number = zeros(1,length(cmap3));
    x_number = zeros(1,length(cmap3));
    s_number = zeros(1,length(cmap3));
    
    for i1 = 1:length(cmap3)
        d = find((rows==i1) + (cols==i1));
        ip_number(i1) = length(find(sum((c(d,:)==1) | (c(d,:)==5),1)));
        p_number(i1) = length(find(sum((c(d,:)==2) | (c(d,:)==6),1)));
        x_number(i1) = length(find(sum(c(d,:)==4,1)));
        s_number(i1) = length(find(sum((c(d,:)==3) | (c(d,:)==7),1)));
    end
    
    dlmwrite(sprintf('./out/p_%s',files(i).name),p_number)
    dlmwrite(sprintf('./out/ip_%s',files(i).name),ip_number)
    dlmwrite(sprintf('./out/x_%s',files(i).name),x_number)
    dlmwrite(sprintf('./out/s_%s',files(i).name),s_number)
    dlmwrite(sprintf('./out/nres_%s',files(i).name),nres)
    dlmwrite(sprintf('./out/coi_%s',files(i).name),coi)
end