function write_convergence(ID_instance, bool)
%
%
%
filename = sprintf('Temp\\SubInstance_%i\\convergence_%i.txt', ID_instance, ID_instance);
fid = fopen(filename, 'w');
fprintf(fid, '%i', bool);
fclose(fid);
end