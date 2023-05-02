function launch_sub_instance(ID_instance, directory)
% 
% /minimize /nosplash /nodesktop if we dont want any GUI
%
subDir = sprintf('Temp\\SubInstance_%i', ID_instance);
scriptfilename = sprintf('%s\\main_sub_instance_%i', subDir, ID_instance);
copyfile('main_sub_instance.m', [scriptfilename, '.m']);
filename = sprintf('%s//instance_%i.bat', subDir, ID_instance);
fid = fopen(filename, 'w');
fprintf(fid, 'cd %s\n', directory);
fprintf(fid, 'cd %s\n', subDir);
fprintf(fid, 'matlab -noFigureWindows -nosplash -nodesktop -r "ID_instance=%i; main_sub_instance_%i"', ID_instance, ID_instance);
fclose(fid);

pause(1)
system(sprintf('Temp\\SubInstance_%i\\instance_%i.bat &', ID_instance, ID_instance));

end