function [new_triangles, new_points, new_N, observer_map, source_map] = QuadBasisFunctionSelect(triangles, basis_select, basis_supports, points)
%Jacques T du Plessis
%September 2019
%19083688@sun.ac.za
% Functionality added for quadrilateral support
%DESCRIPTION


%basis select will recieve the requested observers and sources from the
%user and user and will create a new mesh with only the selected triangle
%required for the new solution to be sent to the mom mex

%The first step in creating the new mesh is selecting and retaining a new
%triangle list that contains only the triangles required for the solution
%in question. This is done by cycling through the basis support list and 
%saving the new triangle list



%variable definitions

new_triangles_selection = [];
new_triangles = [];



basis_select;
new_N = size(basis_select,2);
%make global tri, this is a tri list where the index is the global triangle
%number and the entry is a flag it equals 1 if this triangle is a suppot for
%the selected basis functions and zero if the triangle is not a support

for DOFF_num = 1:size(basis_select,2)

    %go through the triangles support list and populate new triangles with
    %the aforementioned flags
    holder = basis_supports(basis_select(1,DOFF_num),:);
    for select = 1:2
        new_triangles_selection(holder(select)) = 1;
    end
    
end

new_triangles_selection;

for old_index = 1:size(new_triangles_selection,2)

    %now create a new list of triangles that only contain the supports
    %by simply appending the selected triangles to a new list
    if new_triangles_selection(old_index) == 1
        
        new_triangles = [new_triangles;triangles(old_index,:)];
        
    end

end
new_triangles;


%now that we have the remaining triangles we must renumber the triangles as
%wel as the basis functions, but the order must be retained so that the new
%entries in the extracted Zmatrix can be related to the original Zmat
tri = [0, 0, 0];
node_map_selection = [];
for tri_index = 1:size(new_triangles,1)
    
    tri = new_triangles(tri_index,1:4);
    %say which nodes remain
    for iter = 1:4
        node_map_selection(tri(iter)) = 1;
    end
    
end
node_map_selection;


%now the node map must be created for the actual renumbering of the data
new_points = [];
node_map = [];
new_node_index = 0;
for old_node_index = 1:size(node_map_selection,2)

    %if the value in the node map selection is 1 the node apears in the new 
    %mesh and the node must be renumbered
    if node_map_selection(old_node_index) == 1
        new_points = [new_points;points(old_node_index,:)];
        new_node_index = new_node_index + 1;
        node_map(old_node_index) = new_node_index;
        
    end

end
points;
new_points;
node_map;
new_triangles;

%Now use the node map to renumber the triangle nodes
for tri_index = 1:size(new_triangles,1)
    
    %say which nodes remain
    for iter = 1:4
        new_triangles(tri_index,iter) = node_map(new_triangles(tri_index,iter));
    end
    
end 

new_triangles;

%finally we wil renumber the basis functions
%again we create a mapping vector for the basis functions
basis_map_selection = [];
flags = [];
for i = 1:size(basis_select,2)

    basis_map_selection(basis_select(1,i)) = 1;
    flags(1:2,basis_select(1,i)) = basis_select(2:3,i);
    
end
flags;
basis_map_selection;


basis_map = [];
new_basis_index = 0;
for old_basis_index = 1:size(basis_map_selection,2)

    %if the value in the node map selection is 1 the node apears in the new 
    %mesh and the node must be renumbered
    if basis_map_selection(old_basis_index) == 1
        %do three cases
         
            new_basis_index = new_basis_index + 1;
            basis_map(old_basis_index) = new_basis_index;
        
    end

end
basis_map;

observer_map = [];
source_map = [];
num_obs = 0;
num_src = 0;
for i = 1:size(basis_map,2)
    
    
    if basis_map(i) ~= 0 %it is a used basis function
        
        %now we must create the basis function to Zmatrix entry map
        %observers
        if flags(1,i) == 1
           
            %this is an observer and must be mapped
            num_obs = num_obs +1;
            observer_map(basis_map(i)) = num_obs;
            
        elseif flags(1,i) == 0
            observer_map(basis_map(i)) = -1;
        end
        %sources
        if flags(2,i) == 1
           
            
            num_src = num_src +1;
            source_map(basis_map(i)) = num_src;
        elseif flags(2,i) == 0
            
            source_map(basis_map(i)) = -1;
            
        end
    
    end
    
    
end
observer_map;
source_map;
%Now use the basis map to renumber the basis functions
for tri_index = 1:size(new_triangles,1)
    
    
    for iter = 9:12
        
        if (new_triangles(tri_index,iter) ~=-1) && (new_triangles(tri_index,iter) <= size(basis_map,2))
            %renumber the Doffs's
            new_triangles(tri_index,iter) = basis_map(new_triangles(tri_index,iter));
            %take care of the zero
            if new_triangles(tri_index,iter) == 0
                new_triangles(tri_index,iter) = -1;
            end
        elseif  new_triangles(tri_index,iter) >  size(basis_map,2)
            new_triangles(tri_index,iter) = -1;
        end
    end
    
end 
new_triangles;






