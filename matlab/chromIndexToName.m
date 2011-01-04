function name = chromIndexToName(chrom_index)
    name = 'chr';
    if chrom_index <= 18
        name = strcat(name, int2str(chrom_index));
    elseif chrom_index == 19
        name = strcat(name, 'X');
    elseif chrom_index == 20
        name = strcat(name, 'Y');
    end
end

