function test_datafuncs

mesh_fl = load_mesh('circle2000_86_fl');
save_mesh(mesh_fl,'junktest1')
mesh_stnd = load_mesh('circle2000_86_stnd');
save_mesh(mesh_stnd,'junktest2')
mesh_spec = load_mesh('circle2000_86_spec');
save_mesh(mesh_spec,'junktest3')

data = femdata(mesh_fl,0);
save_data(data,'junk.paa');
datasv = load_data('junk.paa');
err = [sum(sum(data.paafl - datasv.paafl)); sum(sum(data.link - datasv.link));...
    sum(sum(data.amplitudefl - datasv.amplitudefl)); sum(sum(data.paax - datasv.paax));...
    sum(sum(data.paaxfl - datasv.paaxfl))]

data = femdata(mesh_stnd,100);
save_data(data,'junk.paa');
datasv = load_data('junk.paa');
err = [sum(sum(data.paa - datasv.paa)); sum(sum(data.link - datasv.link))]

data = femdata(mesh_spec,100);
save_data(data,'junk.paa');
datasv = load_data('junk.paa');
err = [sum(sum(data.paa - datasv.paa)); sum(sum(data.link - datasv.link))]




