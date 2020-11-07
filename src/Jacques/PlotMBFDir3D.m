function [] = PlotMBFDir3D(tri,nodes,node_coords,dof_data,centroid_flowdirs,MBF)

nr_centroids = length(centroid_flowdirs(:,1));
centroid_coords = zeros(nr_centroids ,3);
kk = 0;
for ii = 1:nr_centroids
    for jj = 1:(4-tri)
        kk = kk +1;
        node_coords_new(kk,1) = node_coords(nodes(ii,jj),1);
        node_coords_new(kk,2) = node_coords(nodes(ii,jj),2);
        node_coords_new(kk,3) = node_coords(nodes(ii,jj),3);
        curr_nodes = dof_data.basis_supports(kk,:);
%         if curr_nodes 
        MBF(kk,1) = MBF(ii,jj);
%         centroid_coords(ii,1) = centroid_coords(ii,1) + node_coords(nodes(ii,jj),1); % x
%         centroid_coords(ii,2) = centroid_coords(ii,2) + node_coords(nodes(ii,jj),2); % y
%         centroid_coords(ii,3) = centroid_coords(ii,3) + node_coords(nodes(ii,jj),3); % z
    end
end
PlotMesh(node_coords,nodes, nodes,tri);
% MBF_vals = nonzeros(sum(MBF,2));
% MBF_vals = MBF_vals(1:2:end);
hold on
h = quiver3(node_coords_new(:,1),node_coords_new(:,2),node_coords_new(:,3),zeros(length(MBF),1),zeros(length(MBF),1),MBF(:,1), 2,'LineWidth', 2);

% h = quiver3(centroid_coords(:,1),centroid_coords(:,2),centroid_coords(:,3),centroid_flowdirs(:,1),centroid_flowdirs(:,2),centroid_flowdirs(:,3), 2,'LineWidth', 2);

% If it is necessary to plot the phasor at some other snap in time:
% quiver3(centroid_coords(:,1),centroid_coords(:,2),centroid_coords(:,3),centroid_flowdirs_imag(:,1),centroid_flowdirs_imag(:,2),centroid_flowdirs_imag(:,3), 2,'LineWidth', 2);

end