function PlotMesh(points,element, triangles, tri)

% -------------------------------------------------------------------------
% Quadrilateral plot
% -------------------------------------------------------------------------
node_val = 1;
if tri == 0
    element_points = [points(element(node_val,1),1),points(element(node_val,1),2),points(element(node_val,1),3)];
    % first_val = element_points;
    % Assign coordinate points to each element
    for p = 1:length(element)
        %plot3(points(element(p,1:2),1),points(element(p,1:2),2),points(element(p,1:2),3));
        %plot3(points(element(p,3:4),1),points(element(p,3:4),2),points(element(p,3:4),3));
        % if mod(p,4) == 0
        %     first_val = [points(element(p,1),1),points(element(p,1),2),points(element(p,1),3)];
        % end
        %    if p == ((vertices-1) * 4)
        %        node_val = p - ((vertices-1) * 4);
        %        first_val = [points(element(node_val,1),1),points(element(node_val,1),2),points(element(node_val,1),3)];
        %    end
        temp1 = [points(element(p,1:2),1),points(element(p,1:2),2),points(element(p,1:2),3)];
        temp = [points(element(p,3),1),points(element(p,3),2),points(element(p,3),3)];
        temp2 = [points(element(p,4),1),points(element(p,4),2),points(element(p,4),3)];
        element_points = [element_points;temp1;temp;temp2];
    end

    % Plot 3 sides of the quadrilateral
    % colours = ['b','r','y','g'];
    figure();
    hold on
    for p = 2:4:length(element_points)
        plot3(element_points(p:p+3,1),element_points(p:p+3,2),element_points(p:p+3,3));
    end

    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
end
% -------------------------------------------------------------------------
% Triangle plot
% -------------------------------------------------------------------------

if tri 
    triangle_points = [points(triangles(node_val,1),1),points(triangles(node_val,1),2),points(triangles(node_val,1),3)];
    figure();
    hold on 
    for p = 1:length(triangles)
        temp = [points(triangles(p,1),1),points(triangles(p,1),2),points(triangles(p,1),3)];
        temp1 = [points(triangles(p,1:2),1),points(triangles(p,1:2),2),points(triangles(p,1:2),3)];
        temp2 = [points(triangles(p,3),1),points(triangles(p,3),2),points(triangles(p,3),3)];
        triangle_points = [triangle_points;temp1;temp2;temp];
    end

    % Plot all points on triangle
    hold on
    for p = 2:4:length(triangle_points)-3
       plot3(triangle_points(p:p+3,1),triangle_points(p:p+3,2),triangle_points(p:p+3,3));
    end

    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
end

end