function [] = PlotCurrent3DQuad(InputDataFormat,ShowEdges,quad_nodes,node_coords,quad_currents,centroid_flowdirs)
% Plot the sorface current distribution over a mesh of triangles in 3D. The
% input data (real, scalar-valued current samples) can be provided in
% different formats, which leads to different plotting strategies.
%
% Input arguments:
%
% <InputDataFormat> and <tri_currents>:
%        0  :  <tri_currents> contain the values of the current at the
%              centroids (one value per row)
%        1  :  <tri_currents> contain the values of the current at vertices
%              (three values per row, one row for each triangle, i.e. these
%              are the element-internally evaluated currents, NOT global
%              vertex associated values)  
%
% <tri_nodes> and <node_coords>:
% Mesh description in standard formats.
%
% <centroid_flowdirs>:
% Optional argument specifying the flow direction (3-vector, magnitude
% irrelevant) at the centroids of the triangles.
% 
% 2019-12-15: Created. MMB.

% Init:
maxVal = max(max(quad_currents));
CLim = [0; maxVal];
% CLim = [0; 2.1e-3];
if ShowEdges
    EdgeColor_setting = 'black';
else
    EdgeColor_setting = 'none';
end

% Plot according to the input data format:
figure;
if InputDataFormat == 0
    %trisurf(tri_nodes,node_coords(:,1),node_coords(:,2),node_coords(:,3),tri_currents(:,1));
    patch('Faces',quad_nodes(:,1:4),'Vertices',node_coords,'FaceVertexCData',quad_currents(:,1),'FaceColor','flat','EdgeColor',EdgeColor_setting);  % ,'LineWidth',2
    hold on;
    colormap jet;
    colorbar;
    caxis(CLim);
    axis equal;
elseif InputDataFormat == 1
    % Compose the expandend mesh data, such that current is linearly
    % represented upon each triangle, without interpolation at shared nodes
    % (since in the new represeentation there is no shared nodes)
    numtri           = size(quad_nodes,1);
    numplotnodes     = 4*numtri;
    plotnodecoords   = zeros(numplotnodes,3);
    plotnodecurrents = zeros(numplotnodes,1);
    plotquadnodes     = zeros(numtri,4);
    for ii = 1:numtri
        jj                          = 4*(ii-1) + 1;
        plotnodecoords(jj:jj+3,1:3) = node_coords(quad_nodes(ii,1:4),1:3);
        plotquadnodes(ii,1:4)        = [jj:jj+3];
        plotnodecurrents(jj:jj+3,1) = quad_currents(ii,1:4)';
    end 
    patch('Faces',plotquadnodes,'Vertices',plotnodecoords,'FaceVertexCData',plotnodecurrents,'FaceColor','interp','EdgeColor',EdgeColor_setting);  % ,'LineWidth',2
    hold on;
    colormap jet;
    colorbar;
    caxis(CLim);
    axis equal;
else
    error('Invalid <InputDataFormat>');
end



% set(0,'defaulttextinterpreter','latex')
% % trimesh(TR,,MagnitudeCenterCurrent, );
% 
% %Make a big fancy figure
% figure
% trisurf(TR,MagnitudeCenterCurrent);
% 
% 
% grid off;
% % title('Bowtie Mesh');
% view(-90,90);
% % axis off;