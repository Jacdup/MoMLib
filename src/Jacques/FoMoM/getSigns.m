function [Rho1, Rho2] = getSigns(thisTriangle)

% 'Rho' is a 2x2 matrix containing both the RWG and Linear function magnitudes
% at both vertex points
% [RWG_v1 Lin_v1; RWG_v2 Lin_v2]

sign1 = thisTriangle(1,5);
sign2 =  thisTriangle(1,4); 
sign3 =  thisTriangle(1,10); 
sign4 =  thisTriangle(1,11); 


if sign1*sign2 == 1
%     Rho2 = [-1,-1;-1,1];
    if sign3 < 0
        Rho2 = [-1,1;-1,-1];
    else
       Rho2 = [-1,-1;-1,1]; 
    end
    if sign4 < 0
        Rho1 = [-1,1;-1,-1];
    else
      Rho1 = [-1,-1;-1,1];  
    end
else
    if sign3 < 0
        Rho2 = [1,1;1,-1];
    else
        Rho2 = [1,-1;1,1];
    end
    if sign4 < 0
        Rho1 = [1,1;1,-1];
    else
       Rho1 = [1,-1;1,1]; 
    end
% else
%     Rho1 = [1,1;1,-1];
%     Rho2 = Rho1;
end

end