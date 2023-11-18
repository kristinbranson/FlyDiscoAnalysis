moviefile = '/groups/branson/bransonlab/flydisco_data/VNC2_JRC_SS50996_RigD_20231108T100740/movie.ufmf'

[readfcn,~,~,headerinfo] = get_readframe_fcn(moviefile)
nr = headerinfo.nr
nc = headerinfo.nc
nominal_frame_size = [nr nc] 
nframes = headerinfo.nframes

for i = 31272
    frame = readfcn(i) ;
    if ~isequal(size(frame), nominal_frame_size) ,
        fprintf(2, 'Wrong size for frame %d\n', i) ;
    end
end
