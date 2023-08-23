fname_folder='C:\Users\benja\OneDrive\Documents\MATLAB\h5oina';
fname_slice_data='Mg 1 Specimen 1 SLICE_0156 SLICE_0156 161';

mtex_location='C:\Users\benja\OneDrive\Documents\MATLAB\mtex-5.8.1'; %working with 5.8.1

[ebsd_original,dataset_header,ebsd_patternmatched] = load_h5oina_pm(fname_slice_data,fname_folder,mtex_location);

% optional clean up 
% clear ebsd_o opt_o ebsd1 MSdata2  
% clear fullfile1 fullfile2 location_mtex

%% Plot data
setMTEXpref('xAxisDirection','east');           %aztec
setMTEXpref('zAxisDirection','outofPlane');     %aztec

figure; plot(ebsd_original,ebsd_original.orientations); 
nextAxis
plot(ebsd_patternmatched,ebsd_patternmatched.orientations);