function waiting_status = readWaiting(filename)

fid = fopen(filename, 'r');
waiting_status = fscanf(fid, '%i');
fclose(fid);