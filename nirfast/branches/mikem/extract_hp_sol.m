function sol = extract_hp_sol(solmesh,solution)
% I extract values from hard priors solution meshes!
%
mesh = solmesh;
if exist('solution') == 1
    mesh = read_solution_noplot(solmesh,solution);
end
regs = unique(mesh.region);
for i = 1:length(regs)
    a(i) = sum(mesh.region==regs(i));
    n = find(mesh.region==regs(i));
    node(i) = n(1);
end
[reg_sort ind] = sort(a,'descend');
node = node(ind);
% this matrix is rows of tumor a, tumor b, fg, ad. Columns are as in
% chromscattlist
for i = 1:length(regs)
    sol(i,1:5) = [mesh.conc(node(i),:) mesh.sa(node(i)) mesh.sp(node(i))];
end
end