function [Feko_far_field] = feko_farfield_extract(filename)
% This simply extracts the feko far field from a text file,t uses the text
% you must copy the farfield data from the out file to a text file, see the examples  

fileID = fopen(filename,'r');
% formatSpec = '%f %f %f %f %f %f %*f %*f %*f %*f %*f %*s';
% formatSpec = '%f %f %f %f %f %f %*f %*f %*f %*f %*f %*s';
formatSpec = '%f %f %f %f %f %f %*f %*f %*f %*s';
sizeA = [6 inf];
Feko_far_field = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
Feko_far_field = Feko_far_field';

%size(Z_mat_FEKO)
% Z_mat_abs = abs(Z_mat_FEKO);
% imagesc(Z_mat_abs);
% colorbar;