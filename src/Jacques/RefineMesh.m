function new_points = RefineMesh(Contour, N)
% Description: Creates N new elements for every element in the original contour

 for n = 1:N
     iter = 0;
    for i = 1:length(Contour)-1

        next_point = Contour(i+1,:);
        curr_point = Contour(i,:);
        
        new_point = (next_point + curr_point)/2;
        
        new_points(i+iter:i+1+iter,:) = [curr_point; new_point];
        iter = iter + 1;
    end
    
    new_points(length(Contour)*2 - 1,:) = Contour(end,:);
    Contour = new_points;
 end