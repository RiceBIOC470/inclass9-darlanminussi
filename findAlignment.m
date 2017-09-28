function z = findAlignment(table, blast, N)

table_cell = table2cell(table);

for i = 1:N
    for j = 1:size(blast.Hits,2)
        if table_cell{i,1} == string(blast.Hits(j).Name)
            align = blast.Hits(i).HSPs.Alignment;
            length_align = length(align);
            store_alig = blast.Hits(i).HSPs.Alignment;
            z = ['For ' table_cell{i,1} ':' newline 'Length of alignment is: ' num2str(length_align) newline countAligned(store_alig)];
            disp(z);
        end
    end
end