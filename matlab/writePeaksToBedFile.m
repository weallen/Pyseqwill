function writePeaksToBedFile(fname, peaks)
    fid = fopen(fname, 'w');
    num_lines = length(peaks.chr);
    fmt = '%s\t%d\t%d\tpeak\t999\t+\t%d\t%d\t%s\n';
% chrom chromStart chromEnd name score strand thickStart thickend itemRgb
    for i=1:num_lines
        chr = chromIndexToName(peaks.chr(i));
        start_pos = peaks.start_pos(i);
        end_pos = peaks.end_pos(i);
        state = peaks.state(i);
        color = '255,0,0';
        if state == 1
            color == '255,0,0';
        elseif state == 2
            color = '0,255,0';
        elseif state == 3
            color = '0,0,255';
        end
        fprintf(fid, fmt, chr, start_pos, end_pos, start_pos, end_pos, color);
    end
    fclose(fid);
end
