function dicom2vtk2(fname)
% This function reads a DICOM stack and writes the representation and data
% into a single .vtk file, which could be viewed/manipulated in Paraview.
% fname should be a pattern name that doesn't include the serial numbering
%
% usage: if the series is subject0304_01,subject0304_02,subject0304_03...
%        use dicom2vtk('subject0304');
% outputs: dicoms_fname.vtk
%
% author: venkat krishnawamy/03262010
% last update: 
% part of NIRFAST package
% (C) Hamid Deghani 2008

outfname = ['dicoms_',fname,'.vtk'];
flist=ls([fname '*']);

% extract grid info from the dicom header
dheader = dicominfo(flist(1,:));
w = int32(dheader.Width); % int32 casting is critical here
h = int32(dheader.Height); 
numslices = length(flist);
xres = dheader.PixelSpacing(1);
yres = dheader.PixelSpacing(2);
zres = dheader.SliceThickness;

%read slice data
data = zeros(w*h,(numslices-2)); %pre-allocated data
fprintf('%s','reading slice: '); 
for i = 1:numslices
    temp = dicomread(flist(i,:)); 
    data(:,i) = reshape(temp',w*h,1);    
    fprintf(' %d',i);
end;
data = reshape(data,w*h*numslices,1);

% write VTK file header and an uniform rectilinear grid representation for
% representing the DICOM data
fid = fopen(outfname,'w');

line0 = '# vtk DataFile Version 2.0';
line1 = ['NIRFAST mesh with solutions for case: {', fname,'}'];
line2 = 'ASCII';
line3 = 'DATASET STRUCTURED_POINTS';
fprintf(fid,'%s\n%s\n%s\n',line0,line1,line2,line3);
line4 = ['DIMENSIONS ', num2str(w),' ', num2str(h),' ', num2str(numslices)];
fprintf(fid,'%s\n',line4);
line5 = ['ORIGIN ', ' 0 0 0'];
fprintf(fid,'%s\n',line5);
line6 = ['SPACING ', num2str(xres),' ',num2str(yres),' ',num2str(zres)];
fprintf(fid,'%s\n',line6);

%save dicom data in all slices as scalars defined on the grid nodes
fprintf(fid,'%s\n',['POINT_DATA ',num2str(h*w*numslices)]);
fprintf(fid,'%s\n',['SCALARS ', 'dicom_data', ' unsigned_short 1']);
fprintf(fid,'%s\n','LOOKUP_TABLE default');
fprintf(fid,'%hu\n', data');

fclose(fid);