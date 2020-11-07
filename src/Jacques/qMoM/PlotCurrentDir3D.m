function [] = PlotCurrentDir3D(tri,nodes,node_coords,centroid_flowdirs, centroid_flowdirs_imag)

nr_centroids = length(centroid_flowdirs(:,1));
centroid_coords = zeros(nr_centroids ,3);

for ii = 1:nr_centroids
    for jj = 1:(4-tri)
        centroid_coords(ii,1) = centroid_coords(ii,1) + node_coords(nodes(ii,jj),1); % x
        centroid_coords(ii,2) = centroid_coords(ii,2) + node_coords(nodes(ii,jj),2); % y
        centroid_coords(ii,3) = centroid_coords(ii,3) + node_coords(nodes(ii,jj),3); % z
    end
    if tri
        centroid_coords(ii,:) = (1/3)*(centroid_coords(ii,:));
    else
        centroid_coords(ii,:) = 0.25*(centroid_coords(ii,:));
    end
end
% PlotMesh(node_coords,nodes, nodes,tri);
hold on
h = quiver3(centroid_coords(:,1),centroid_coords(:,2),centroid_coords(:,3),centroid_flowdirs(:,1),centroid_flowdirs(:,2),centroid_flowdirs(:,3), 2,'LineWidth', 2);

% If it is necessary to plot the phasor at some other snap in time:
% quiver3(centroid_coords(:,1),centroid_coords(:,2),centroid_coords(:,3),centroid_flowdirs_imag(:,1),centroid_flowdirs_imag(:,2),centroid_flowdirs_imag(:,3), 2,'LineWidth', 2);

end