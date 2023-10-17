function [EBSDpat]=loadEBSP_T(fname,dataset_header,pattern_number)
%% 


%% 
EBSP_PW=double(dataset_header.Pattern_Width);
EBSP_PH=double(dataset_header.Pattern_Height);
% loc=header.fileloc;
% fname=header.fname;
Pat = [h5info(fname).Groups(1).Name '/EBSD/Data/Processed Patterns'];

EBSDpat=double((rot90(shiftdim(h5read(fname,Pat,[1 1 pattern_number],[EBSP_PW EBSP_PH 1])))));

%have to install the filter, as per the asnwer here
% https://www.mathworks.com/matlabcentral/answers/393621-how-can-i-import-a-lzf-compressed-hdf5-dataset-in-matlab
%







