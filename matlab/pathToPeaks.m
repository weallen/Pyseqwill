function peaks = pathToPeaks(path, windows)
    num_state_trans = 0;
    curr_state = path(1);
    for i=2:length(path)
        if path(i) ~= curr_state
            num_state_trans = num_state_trans + 1;
            curr_state = path(i);
        end
    end
    peaks = struct('chr', zeros(1, num_state_trans), ...
                    'start_pos', zeros(1, num_state_trans), ...
                    'end_pos', zeros(1, num_state_trans), ...
                    'state', zeros(1, num_state_trans));
    prev_state = path(1);
    prev_start = 1;
    peak_idx = 1;
    for i=2:length(path)
        curr_state = path(i);
        if prev_state ~= curr_state || windows.chr(prev_start) ~= windows.chr(i)
            % split into 2 if crosses chromosome boundary
            peaks.chr(peak_idx) = windows.chr(prev_start);
            peaks.start_pos(peak_idx) = windows.start(prev_start);
            peaks.end_pos(peak_idx) = windows.xEnd(i-1);
            peaks.state(peak_idx) = curr_state;
            prev_start = i;
            peak_idx = peak_idx + 1;
        end    
    end
end
