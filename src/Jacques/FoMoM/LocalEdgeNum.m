function [local_edge_num] = LocalEdgeNum(edge_nodes,tri_nodes)
% The local edge number is returned. Is is assumed that the input arguments
% are sorted lists of global node numbers, of lengths 2 and 3, respectively.
%
% 2019-12-13: Created. MMB.

%edge_select        = [1 2 ;1 3 ;2 3];
%direction_select   = [6 5 4];
%doff_select        = [9 8 7];

% edge_nodes      = [1 2]    = curr_edge;
% tri_nodes(t1,:) = [1 2 27] = global triangle nodes


local_edge_nodes_def = [1 2
    1 3
    2 3];


% if tri_nodes(2) > tri_nodes(3)% Then it's not sorted, and numbering should change
%     if edge_nodes == tri_nodes(local_edge_nodes_def(1,:))
%         %     local_edge_num = 1;
%         local_edge_num = [2 1];
%         return;
%     elseif edge_nodes == tri_nodes(local_edge_nodes_def(2,:))
%         %     local_edge_num = 2;
%         local_edge_num = [1 2];
%         return;
%     else
%         local_edge_num = [2 1];
%         %     local_edge_num = 1;
%     end
%     
% else
    
    if edge_nodes == tri_nodes(local_edge_nodes_def(1,:))
        %     local_edge_num = 1;
        local_edge_num = [1 2];
        return;
    elseif edge_nodes == tri_nodes(local_edge_nodes_def(2,:))
        %     local_edge_num = 2;
        local_edge_num = [2 1];
        return;
    else
        local_edge_num = [1 2];
        %     local_edge_num = 1;
    end
% end