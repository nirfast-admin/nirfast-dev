function batch_test_newlink

mesh_fl = 'circle2000_86_fl';
mesh_stnd = 'circle2000_86_stnd';
mesh_spec = 'circle2000_86_spec';
err = [];
err = [err; test_new_link_spec(mesh_spec)];
err = [err; test_new_link_fl(mesh_fl)];
err = [err; test_new_link_stnd(mesh_stnd)];