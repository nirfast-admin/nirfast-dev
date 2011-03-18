function [data] = combine_cursor_data(data,mesh,index)
% Makes a cell of combined cursor data for removing and then removes all
% marked data points. Must be selected and removed using data cursor in
% each figure.


%%%%%%%%%%%%%%%%%%
%index = [3 4 5 6]; % set me based on which wavelengths were dropped
%%%%%%%%%%%%%%%%%%
cursor = {0 0 0 0 0 0};
for i = 1:6
    if any(index == i)
       cursor{i} = evalin('base',['c',num2str(i)]); 
    end
end
for i = index
    [data] = remove_cursor_data_func(cursor{i},data,mesh,i);
end
clear index cursor
