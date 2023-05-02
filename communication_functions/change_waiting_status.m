function change_waiting_status(ID_instance)
filename = sprintf('Temp\\SubInstance_%i\\%s_%i.txt', ID_instance, 'waiting_file', ID_instance);

fid = fopen(filename, 'r');
waiting_status = fscanf(fid, '%i');
fclose(fid);

fid = fopen(filename, 'w');
fprintf(fid, '%i', ~waiting_status);
fclose(fid);
end