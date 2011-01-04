function writePathToWigFile(fname, path, windows, span)
    fid = fopen(fname, 'w');
    num_lines = length(path(path == 2));
    chr_idx = zeros(1,21);
    for i=1:20
        chr_idx(i) = min(find(windows.chr == i));
    end
    chr_idx(21) = length(windows.chr)+1;
    fprintf(fid, 'track type=wiggle_0\n');
    for i=1:20
        chr = chromIndexToName(i);
        fprintf(fid, 'variableStep chrom=%s span=%d\n', chr, span);
        chr_start = chr_idx(i);
        chr_stop = chr_idx(i+1)-1;
        chr_len = chr_stop - chr_start;
        for j=1:chr_len
            if path(chr_start + j) == 2
                fprintf(fid, '%d 10\n', windows.start(chr_start + j));
            end
        end
    end
    fclose(fid);
end

