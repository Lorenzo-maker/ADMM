function convergence = read_convergence(ID_instance)
%
%
%
filename = sprintf('SubInstance_%i\\convergence_%i.txt', ID_instance, ID_instance);
fid = fopen(filename, 'r');
convergence = fscanf(fid, '%i');
fclose(fid);
end