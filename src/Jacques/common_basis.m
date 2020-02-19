function [all_basis_functions, num_obs, num_src] = common_basis(obs_basis_select_cells,src_basis_select_cells)
%Robey Beswick
%Februaury 2019
%18472648@sun.ac.za
%
%DESCRIPTION
%This function returns a 3 by N vector of the union of the basis functions 
%selected, it is recieved by basis function select and is used for three
%functions, first we use this to find the basis functions that are relevant
%to the new mesh, second for reumbering the basis functions in order and
%third for creating maps that allow us to map the basis functions to the
%new entries in the Zmatrix.

all_basis_functions = [];
obs_basis_select = [];
src_basis_select = [];

%convert the cell matrix to number matix
for index = 1:size(obs_basis_select_cells,1)
    
    obs_basis_select = [obs_basis_select,cell2mat(obs_basis_select_cells(index))];
end
num_obs = size(obs_basis_select,2);
%convert the cell matrix to number matix
for index = 1:size(src_basis_select_cells,1)
    
    src_basis_select = [src_basis_select,cell2mat(src_basis_select_cells(index))];
end
num_src = size(src_basis_select,2);

%Now we must ind the common entries to the observers and cells
i = 1;
j = 1;

%Note that this function does not handle repeated entries
while (i <= size(obs_basis_select,2)) && (j <= size(src_basis_select,2))
    
   if obs_basis_select(i) <  src_basis_select(j)
       %value is only in observer save it an go to next entry
       hold = [obs_basis_select(i);1;0];
       all_basis_functions = [all_basis_functions,hold];
       i = i + 1;
       
   elseif src_basis_select(j) < obs_basis_select(i)
       %value is only in array 2 save it and go to next entry
       hold = [src_basis_select(j);0;1];
       all_basis_functions = [all_basis_functions,hold];
       j = j + 1;
   else
       hold = [obs_basis_select(i);1;1];
       all_basis_functions = [all_basis_functions,hold];
       i = i + 1;
       j = j + 1;
   end
end

%The elements remaining in the larger array must must still be printed
while i <= size(obs_basis_select,2)
    hold = [obs_basis_select(i);1;0];
       all_basis_functions = [all_basis_functions,hold];
       i = i + 1;
end

while j <= size(src_basis_select,2)
    hold = [src_basis_select(j);0;1];
       all_basis_functions = [all_basis_functions,hold];
       j = j + 1;
end

