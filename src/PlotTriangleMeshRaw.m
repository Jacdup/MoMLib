function [] = PlotTriangleMeshRaw(node_coords,tri_nodes,TextOn)
% Plot a mesh of triangles with node and triangle numbers.
% node_coords    : M x (X,Y,Z)
% tri_nodes      : N x (node1,node2,node3)
%
% 2019-12-11: Created. MMB. 

% Init:
numnodes = size(node_coords,1);
numtri   = size(tri_nodes,1);

% Plot the nodes:
figure;
plot3(node_coords(:,1),node_coords(:,2),node_coords(:,3),'ko');
hold on;
axis equal;
%axis off;
if TextOn 
    for ii = 1:numnodes
        text(node_coords(ii,1),node_coords(ii,2),node_coords(ii,3),[' ',num2str(ii)],'HorizontalAlignment','left','Color','k'); 
    end
end

% Plot the triangles:
%trimesh(tri_nodes,node_coords(:,1),node_coords(:,2),node_coords(:,3),'EdgeColor','blue','FaceColor','none');
patch('Faces',tri_nodes,'Vertices',node_coords,'EdgeColor','blue','FaceColor','none');  % ,'LineWidth',2
% for ii = 1:numtri
%     n1 =  tri_nodes(ii,1);
%     n2 =  tri_nodes(ii,2);
%     n3 =  tri_nodes(ii,3);
%     plot3([ node_coords(n1,1) node_coords(n2,1) ],[ node_coords(n1,2) node_coords(n2,2) ],[ node_coords(n1,3) node_coords(n2,3) ],'b');
%     plot3([ node_coords(n1,1) node_coords(n3,1) ],[ node_coords(n1,2) node_coords(n3,2) ],[ node_coords(n1,3) node_coords(n3,3) ],'b');
%     plot3([ node_coords(n2,1) node_coords(n3,1) ],[ node_coords(n2,2) node_coords(n3,2) ],[ node_coords(n2,3) node_coords(n3,3) ],'b');
% end
if TextOn 
    for ii = 1:numtri
        n1 =  tri_nodes(ii,1);
        n2 =  tri_nodes(ii,2);
        n3 =  tri_nodes(ii,3);
        tricentroid = (1/3) * (node_coords(n1,:) + node_coords(n2,:) + node_coords(n3,:));
        text(tricentroid(1),tricentroid(2),tricentroid(3),num2str(ii),'HorizontalAlignment','center','Color','b'); % ['t',num2str(ii)]
    end
end


%scatter3(array_points(:,1), array_points(:,2), array_points(:,3));
% figure()
% PlotTriangles(array_points, array_triangles(:,1:3))


%function [] = PlotTriangles(points, triangles)
%TR = triangulation(triangles,points);
%trimesh(TR);