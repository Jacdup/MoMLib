function [] = PlotTriangleMeshDofs(mesh_data,TextOn)
% Plot a mesh of triangles with node and triangle and global edge numbers.
% mesh_data    : output of the function <CreateMeshData.m>
%
% 2019-12-12: Created. MMB. 

% Plot the nodes and triangles:
PlotTriangleMeshRaw(mesh_data.node_coords,mesh_data.tri_nodes,TextOn);

% Get edge midpoints:
numedges       = size(mesh_data.edges,1);
ALPHA          = 0.5;
edge_midpoints = ALPHA*mesh_data.node_coords(mesh_data.edges(:,1),:) + (1-ALPHA)*mesh_data.node_coords(mesh_data.edges(:,2),:);
plot3(edge_midpoints(:,1),edge_midpoints(:,2),edge_midpoints(:,3),'m*');
if TextOn 
    for ii = 1:numedges
        text(edge_midpoints(ii,1),edge_midpoints(ii,2),edge_midpoints(ii,3),[' ',num2str(ii)],'HorizontalAlignment','left','Color','m');
    end
end


% % Plot the global edges:
% for ii = 1:numedges
%     thisedge =  mesh_data.edges(ii,1:2);
%     %plot3( mesh_data.node_coords(thisedge,1), mesh_data.node_coords(thisedge,2), mesh_data.node_coords(thisedge,3), 'm--');
%     ALPHA = 0.25;
%     tempcoords = ALPHA*mesh_data.node_coords(thisedge(1),:) + (1-ALPHA)*mesh_data.node_coords(thisedge(2),:);
%     plot3(tempcoords(1), tempcoords(2), tempcoords(3), 'm*');
%     tempcoords = ALPHA*mesh_data.node_coords(thisedge(2),:) + (1-ALPHA)*mesh_data.node_coords(thisedge(1),:);
%     plot3(tempcoords(1), tempcoords(2), tempcoords(3), 'm*');
%     if TextOn 
%         edgecentroid = (1/2) * (mesh_data.node_coords(thisedge(1),:) + mesh_data.node_coords(thisedge(2),:));
%         text(edgecentroid(1),edgecentroid(2),edgecentroid(3),num2str(ii),'HorizontalAlignment','center','Color','m');
%     end
% end
