function [data] = remove_cursor_data(cursor,data,wv_index)
% removes cursor data for a 6wv data set.
% assumes that cursor is a struct with fields c1, c2, c3... Made from
% combine_cursor_data.m

data = remove_cursor_data_func(cursor.c1,data,1);
data = remove_cursor_data_func(cursor.c2,data,2);
data = remove_cursor_data_func(cursor.c3,data,3);
data = remove_cursor_data_func(cursor.c4,data,4);
data = remove_cursor_data_func(cursor.c5,data,5);
data = remove_cursor_data_func(cursor.c6,data,6);
end    
    
    