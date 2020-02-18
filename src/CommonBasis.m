function [all_basis_functions, num_obs, num_src, global_to_redu_dofs] = CommonBasis(obs_basis_select,src_basis_select,num_dofs,SortCheckOn)
% This function returns a 3 by N vector of the union of the basis functions 
% selected, it is recieved by basis function select and is used for three
% functions, first we use this to find the basis functions that are relevant
% to the new mesh, second for renumbering the basis functions in order and
% third for creating maps that allow us to map the basis functions to the
% new entries in the Zmatrix.
%
% Structure of a column (num col = num total common basis functions --- union of the obs and src lists)
% [  full_matrix_dof_number   observer_flag   source_flag  ]'
% 
%
% Structure of <global_to_redu_dofs>: number of cols is equal to number of
% original ("global") dofs. Entries are either -1, or the corresponding dof
% number in the reduced selection (i.e. column number in <all_basis_functions>).
%
% This routine only accepts sorted input observer and source indices
% (although this check can be turned off for speed). It is only
% suitable for calculating a collapsed block of the MoM 
% matrix. In this collapsed block (e.g. obs [10 51] and src [2 3 235] will
% eventually yield a 2 x 3 matrix), increasing collapsed indices maps to
% increasing global indices. Therefore, this routine is not an avenue for
% independent obs and src permutations, which will generally yield
% unsymmetric diagonal blocks (if that is what was selected). A future
% extension could be to support that.
%
% 2019-02-00: Created. Robey Beswick, 18472648@sun.ac.za.
% 2019-12-13: Polished up, added error check, removed unncecesary use of
% cell array input arguments. MMB

% Initialise the outputs:
all_basis_functions = [];
num_obs             = size(obs_basis_select,2);
num_src             = size(src_basis_select,2);
global_to_redu_dofs = -ones(1,num_dofs);

% Check that the vectors are sorted in ascending order:
if SortCheckOn 
    if obs_basis_select ~= sort(obs_basis_select,2) 
        error('Invalid obs order in CommonBasis');
    end
    if src_basis_select ~= sort(src_basis_select,2) 
        error('Invalid src order in CommonBasis');
    end
end

%Now we must find the common entries to the observers and cells:
%(Note that this function does not handle repeated entries)
i = 1;
j = 1;
while (i <= num_obs) && (j <= num_src)
    
   if obs_basis_select(i) < src_basis_select(j)
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

%The elements remaining in the larger array must still be printed
while i <= num_obs
    hold = [obs_basis_select(i);1;0];
       all_basis_functions = [all_basis_functions,hold];
       i = i + 1;
end

while j <= num_src
    hold = [src_basis_select(j);0;1];
       all_basis_functions = [all_basis_functions,hold];
       j = j + 1;
end

% Fill the forward map:
for ii = 1:size(all_basis_functions,2)
    global_to_redu_dofs(1,all_basis_functions(1,ii)) = ii;
end

disp('<observer_map> and <source_map> should be constructed in CommonBasis, rather than in BasisFunctionSelect');
