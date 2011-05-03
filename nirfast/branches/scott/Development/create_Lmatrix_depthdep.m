function create_Lmatrix_depthdep(Jsf,nnodes,uniqreg,meshregion)


tic
L = sparse(eye(nnodes))*Jsf;

for i = 1:length(uniqreg)
    R = find(meshregion == uniqreg(i));
    ty = -1/(length(R)+1);
    L(R,R) = ty;
end

d = -diag(L)+1;
L = (L + diag(d));
L_paa = sparse([L zeros(size(L)); zeros(size(L)) L]);

time = toc;
disp(['Computation Time L-matrix = ' num2str(time) ]); 