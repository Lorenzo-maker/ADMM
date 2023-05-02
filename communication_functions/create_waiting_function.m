function filename = create_waiting_function(ID_instance)

filename = sprintf('Temp\\SubInstance_%i\\waiting_file_%i.txt', ID_instance, ID_instance);
fid = fopen(filename, 'w');
fprintf(fid, '0');
fclose(fid);


end