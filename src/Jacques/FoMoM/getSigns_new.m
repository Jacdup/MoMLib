function [Rho, Rho2] = getSigns_new(triangle,col1,col2, i,i2)
% 
% if i == 680
%     test =1;
% end

% if i2 == 37
%     test = 1;
% end
% 
% if i == i2
%     test =1;
% end

Rho(1,1) = triangle(i,col1-3);
Rho(1,2) = 0.5* triangle(i,col1+3);
Rho(2,1) = triangle(i,col1-3);
Rho(2,2) = Rho(1,2) * -1;
% Rho(2,2) = 0.5* triangle(i,col1+3);
% Rho(1,2) =  Rho(2,2) * -1;

Rho2(1,1) = triangle(i2,col2-3);
Rho2(1,2) = 0.5* triangle(i2,col2+3);
Rho2(2,1) = triangle(i2,col2-3);
Rho2(2,2) = Rho2(1,2) * -1;

% Rho2(2,2) =  0.5* triangle(i2,col2+3);
% Rho2(1,2) = Rho2(2,2) * -1;

% Rho(1,1) = Rho(1,1) * -1;
% Rho2(1,1) = Rho2(1,1) * -1;
% Rho(2,1) = Rho(2,1) * -1;
% Rho2(2,1) = Rho2(2,1) * -1;

% Rho(1,1) = 1; % Do nothing with RWG
% Rho2(1,1) = 1;
% Rho(2,1) = 1;
% Rho2(2,1) = 1;
            
            
%             Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col+extra_dof_col);
%             Rho(1,1,Rho_index) = Rho(1,1,Rho_index) * triangle_blah(i-1,4);
%             Rho(1,2,Rho_index) = Rho(1,2,Rho_index) * 0.5* triangle_blah(i-1,10);
%             Rho(2,1,Rho_index) = Rho(2,1,Rho_index) * triangle_blah(i-1,4);
%             Rho(2,2,Rho_index) = Rho(1,2,Rho_index) * -1;