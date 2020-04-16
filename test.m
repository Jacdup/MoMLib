
for MBF_num = 1:3
    col = 1;
%     col_index = MBF_num-1;
    for MBF_node = 1:3
        
        col_index = col + (MBF_num-1)
        col = col + 3;
%         col_index = ((MBF_num-1)*(3)) + (MBF_node)
%         col_index = MBF_num+MBF_node
    end
end