% 
% 
% median_multi_rep_spectra('Basis_903nm_trans','Median_data_903nm.txt',10,0)
% basis = create_basis_spectra('Median_data_903nm.txt')
% basis(find(basis<0))=0;
% save('Basis_903nm.txt','basis','-ASCII');

median_multi_rep_spectra('Basis_912nm_trans','Median_data_912nm',10,0)
basis = create_basis_spectra('Median_data_912nm.raw');
basis(find(basis<0))=0;
save('Basis_912nm.txt','basis','-ASCII');
