function [points, triangles] = ImportTriangleMeshNastran(filename)%filename recieved is a string
% Import a mesh of triangles from a Nastran *.nas file. The tested Nastran files
% came from FEKO, where the FEKO mesh was exported to a *.nas file.
%
% 2018-09-00: Created. Robey Beswick 18472648@sun.ac.za.
% 2019-12-11: MMB renamed and 

%clear all;

fileID = fopen(filename,'r');
tline = 0;
points = [];
triangles = [];

while tline ~= -1
    tline = fgetl(fileID);
    %first fetch the co-ordinate points
    if tline(1) == 'G'
        co_ord_str = tline;
    elseif tline(1) == '*'%Third coordinate from new line
        co_ord_str = [co_ord_str,tline];
        %Now remove the numbers from the string and convert them to floats
        points = [points;sscanf(co_ord_str,'%*7c %*d %f %f %*d %*4c %*d %f')'];
        co_ord_str = 0;      
    end
    %Second we fetch the node data
    if tline(1) == 'C'
        next_tri  = sscanf(tline,'%*7c %*d %*d %d %d %d')';
        triangles = [triangles;next_tri];
    end
end
fclose(fileID);
%triangles = cell2mat(triangles);  % this is just simply because importnastran had returned a list of cells previously



