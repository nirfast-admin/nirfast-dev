function [paa] = dataprep_spec(fn,rep)

% read spectrometer data.


tic
close all;
%fn = 'inclusion_gel'; % base before any labview added things.
%rep = 1;

% load log file
[log1] = load_log_file([fn,'_side1_trans']);
[log2] = load_log_file([fn,'_side2_trans']);

for i = find(log1.sources==1)'
    data{i} = importdata([fn,'_side1_trans_s',num2str(i),'_rep',num2str(rep),'.raw']);
end
for i = find(log2.sources==1)'
    data{i} = importdata([fn,'_side2_trans_s',num2str(i),'_rep',num2str(rep),'.raw']);
end

% find median from pixels 40-200 and subtract from all data points.
for i = [find(log1.sources==1)' find(log2.sources==1)']
    med = median(data{i}(40:220,:));
    med = repmat(med,1340,1);
    data{i} = data{i} - med;
    %plot(data{i}(:,10))
end


% load basis spectra
basis903 = load('Basis_903nm.txt');
basis912 = load('Basis_912nm3.txt');

% smooth and normalize
for i = 1:16
    basis903(:,i) = smooth(basis903(:,i));
    basis903(:,i) = basis903(:,i)/max(basis903(:,i));
    basis912(:,i) = smooth(basis912(:,i));
    basis912(:,i) = basis912(:,i)/max(basis912(:,i));
end

% fit side 1 and integrate
ii = 1;
for i = find(log1.sources==1)'
    for j = find(log1.meas==1)'
        close all
        clear error
        % smooth data using median filter
        data{i}(:,j) = smooth(data{i}(:,j),10);
        data_max  = max(data{i}(:,j));
        % find peaks
        [peakval,loc] = findpeaks(data{i}(:,j),'minpeakheight',data_max/2);
        %plot(data{i}(:,j))
        %hold on
        %plot(loc,data{i}(loc,j),'ro')
        b = find(abs(loc(end)-loc)<20);
        a = find(abs(loc(1)-loc)<20);
        b = round(mean(loc(b)));
        a = round(mean(loc(a)));
        % pick range +/- 50 from peaks for fitting
        aa = a-50; bb = b+50;
        
        % scan basis 912 for best fit
        basis912_max = round(median(find(basis912(:,j)==1)));
        basis912(:,j) = circshift(basis912(:,j),b-basis912_max);
        clear error
        for k = 1:7
            basis2 = circshift(basis912(:,j),k-4);
            [spec error(k)] = LLS_specfit([basis903(aa:bb,j) basis2(aa:bb)],data{i}(aa:bb,j),0);
        end
        [mina minb] = min(error);
        basis2 = circshift(basis912(:,j),minb-4);
        
        % scan basis 903 for best fit
        basis903_max = round(median(find(basis903(:,j)==1)));
        basis903(:,j) = circshift(basis903(:,j),a-basis903_max);
        clear error
        for k = 1:7
            basis1 = circshift(basis903(:,j),k-4);
            [spec error(k)] = LLS_specfit([basis1(aa:bb) basis2(aa:bb)],data{i}(aa:bb,j),0);
        end
        [mina minb] = min(error);
        basis1 = circshift(basis903(:,j),minb-4);
        [spec e] = LLS_specfit([basis1(aa:bb) basis2(aa:bb)],data{i}(aa:bb,j),1);
        
        % find peak 948
        data948 = data{i}(b+100:b+300,j);
        med = median(data948(1:30));
        [peakval,loc] = findpeaks(data948,'minpeakheight',med*5); % SNR=5
        loc = round(median(loc));
%         plot(data948)
%         hold on
%         plot(loc,data948(loc),'ro')
        

        % integrate each peak and build data file
        side1(ii,1:2) = trapz(spec);
        side1(ii,3) = trapz(data948(loc-40:loc+40));
        ii = ii+1;
    end
end


% fit side 2 and integrate
ii = 1;
for i = find(log2.sources==1)'
    for j = find(log2.meas==1)'
        close all
        clear error
        % smooth data using median filter
        data{i}(:,j) = smooth(data{i}(:,j),10);
        data_max  = max(data{i}(:,j));
        % find peaks
        [peakval,loc] = findpeaks(data{i}(:,j),'minpeakheight',data_max/2);
         plot(data{i}(:,j))
        hold on
         plot(loc,data{i}(loc,j),'ro')
        b = find(abs(loc(end)-loc)<20);
        a = find(abs(loc(1)-loc)<20);
        b = round(mean(loc(b)));
        a = round(mean(loc(a)));
        % pick range +/- 50 from peaks for fitting
        aa = a-50; bb = b+50;
        
        % scan basis 912 for best fit
        basis912_max = round(median(find(basis912(:,j)==1)));
        basis912(:,j) = circshift(basis912(:,j),b-basis912_max);
        clear error
        for k = 1:7
            basis2 = circshift(basis912(:,j),k-4);
            [spec error(k)] = LLS_specfit([basis903(aa:bb,j) basis2(aa:bb)],data{i}(aa:bb,j),0);
        end
        [mina minb] = min(error);
        basis2 = circshift(basis912(:,j),minb-4);
        
        % scan basis 903 for best fit
        basis903_max = round(median(find(basis903(:,j)==1)));
        basis903(:,j) = circshift(basis903(:,j),a-basis903_max);
        clear error
        for k = 1:7
            basis1 = circshift(basis903(:,j),k-4);
            [spec error(k)] = LLS_specfit([basis1(aa:bb) basis2(aa:bb)],data{i}(aa:bb,j),0);
        end
        [mina minb] = min(error);
        basis1 = circshift(basis903(:,j),minb-4);
        [spec e] = LLS_specfit([basis1(aa:bb) basis2(aa:bb)],data{i}(aa:bb,j),0);
        
        % find peak 948
        data948 = data{i}(b+100:b+300,j);
        med = median(data948(1:30));
        [peakval,loc] = findpeaks(data948,'minpeakheight',med*5); % SNR=5
        loc = round(median(loc));
%         plot(data948)
%         hold on
%         plot(loc,data948(loc),'ro')
        

        % integrate each peak and build data file
        side2(ii,1:2) = trapz(spec);
        if ~isnan(loc)
            side2(ii,3) = trapz(data948(loc-40:loc+40));
        end
        ii = ii+1;
    end
end


% build data file using n sources, m detectors.
paa.wv = [903;912;948];
data_tmp = [side1; side2];

% side 1
k = 1;
kk = 1;
for i = 1:8
    for j = 1:16
        paa.link(k,1) = i; paa.link(k,2) = j;
        if log1.exp_times(j,i)~=0
            paa.link(k,3:5) = 1;
            paa.paa(k,1:3) = data_tmp(kk,:)./log1.exp_times(j,i);
            k=k+1; kk=kk+1;
        else
            paa.link(k,3:5) = 0;
            paa.paa(k,1:3) = 0;
            k=k+1;
        end
    end
end
% side 2
for i = 9:16
    for j = 1:16
        paa.link(k,1) = i; paa.link(k,2) = j;
        if log2.exp_times(j,i)~=0
            paa.link(k,3:5) = 1;
            paa.paa(k,1:3) = data_tmp(kk,:)./log2.exp_times(j,i);
            k=k+1; kk=kk+1;
        else
            paa.link(k,3:5) = 0;
            paa.paa(k,1:3) = 0;
            k=k+1;
        end
    end
end

% for safe keeping
% % side 1
% k = 1;
% kk = 1;
% for i = 1:8
%     for j = 1:16
%         paa.link(k,1) = i; paa.link(k,2) = j;
%         if (log1.sources(i)==1 && log1.meas(j)==1)
%             paa.link(k,3:5) = 1;
%             paa.paa(k,1:3) = data_tmp(kk,:)./log1.exp_times(j,i);
%             k=k+1; kk=kk+1;
%         else
%             paa.link(k,3:5) = 0;
%             paa.paa(k,1:3) = 0;
%             k=k+1;
%         end
%     end
% end
% % side 2
% for i = 9:16
%     for j = 1:16
%         paa.link(k,1) = i; paa.link(k,2) = j;
%         if (log2.sources(i)==1 && log2.meas(j)==1)
%             paa.link(k,3:5) = 1;
%             paa.paa(k,1:3) = data_tmp(kk,:)./log2.exp_times(j,i);
%             k=k+1; kk=kk+1;
%         else
%             paa.link(k,3:5) = 0;
%             paa.paa(k,1:3) = 0;
%             k=k+1;
%         end
%     end
% end
close all
toc





