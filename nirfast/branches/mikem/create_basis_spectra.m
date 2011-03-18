function basis = create_basis_spectra(fn_median_data)

param = eval('input_calib');
basis = load(fn_median_data);

for i = 1:16
    if i ~= 3 && i ~=9
        zero_range_data = basis(param.fl_baseline_region_pixels,i);
        basis(:,i) = basis(:,i)-median(zero_range_data);
        basis(:,i) = basis(:,i)/max(basis(:,i));
    end
end
