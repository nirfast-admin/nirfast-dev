function [spec error] = LLS_specfit(A,data,flag)

% Performs least squares fitting between measured spectrum (data), and
% measured basis spectra (in matrix A).

% input basis spectra in matrix A with columns denoting spectra.
% input data as column vector
% output data is matrix A with each column scaled by its minimization
% coefficient

[npix,m] = size(A);

H = A'*A;
b = A'*data;
u = H\b;

% Constrain negative fits
% for i = 1:m
%     if u(i) < 0
%         act(i) = 0;
%     else
%         act(i) = 1;
%     end
% end
% 
% if sum(act) < m
%     ind = find(act == 0);
%     A2 = A;
%     A2(:,ind) = [];
%     H = A2'*A2;
%     b = A2'*data;
%     u = H\b;
%     k = 1;
%     for i = 1:m
%         if act(i) == 1
%             u2(i) = u(k);
%             k = k+1;
%         else
%             u2(npix,i) = 0;
%         end
%     end
%     u = u2;
% end

% Collect unmixed spectra
spec=[];
for i = 1:m
    spec = [spec, u(i)*A(:,i)];
end
clear error
% compute error
fit = spec(:,1)+spec(:,2);
error = mean((data - fit).^2);


% Plot flag
if flag == 1
    figure
    x_axis = [1:npix]';
    plot(x_axis,data, '-k', x_axis,spec(:,1),x_axis,spec(:,2),x_axis,fit,':');
    ylabel('counts/s'); xlabel('pixel number');
end

