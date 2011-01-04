function data = loadWindows( h5filename )
    loc = H5F.open(h5filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    try
        data = loadWindowsData(loc, 'windows');
        H5F.close(loc);
    catch exc
        H5F.close(loc);
        rethrow(exc);
    end
end

%class GenomeWindow:
%    idx = tables.Int32Col() # absolute index of window
%    chr = tables.Int8Col() # chromosome 
%    start = tables.Int32Col() # end of window
%    end = tables.Int32Col() # start of window
%    chr_pos = tables.Int32Col()

function data = loadWindowsData(loc, subpath)
    int32_type = H5T.copy('H5T_NATIVE_INT32');
    sz(1) = H5T.get_size(int32_type);
    int8_type = H5T.copy('H5T_NATIVE_INT8');
    sz(2) = H5T.get_size(int8_type);
    sz(3) = H5T.get_size(int32_type);
    sz(4) = H5T.get_size(int32_type);
    sz(5) = H5T.get_size(int32_type);
    offset(1) = 0;
    offset(2:5) = cumsum(sz(1:4));
    dset = H5D.open(loc, subpath); 
    win_type = H5T.create('H5T_COMPOUND', sum(sz));
    H5T.insert(win_type, 'idx', offset(1), int32_type); 
    H5T.insert(win_type, 'chr', offset(2), int8_type);
    H5T.insert(win_type, 'start', offset(3), int32_type);
    H5T.insert(win_type, 'end', offset(4), int32_type);
    H5T.insert(win_type, 'chr_pos', offset(5), int32_type);
    space = H5D.get_space(dset);
    [~, dims, ~] = H5S.get_simple_extent_dims (space);
    dims = fliplr(dims);
    data = H5D.read(dset, win_type, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
    H5D.close(dset);
end
