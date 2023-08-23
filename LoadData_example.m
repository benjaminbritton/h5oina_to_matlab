% Example use of the h5oinA loading script

%folder that contains the h5oina files, 
% e.g. fname_folder='C:\Users\benja\OneDrive\Documents\MATLAB\h5oina';
fname_folder='C:\Users\benja\OneDrive\Documents\MATLAB\h5oina';

%h5oina file, excluding final file extension & folder
%if the '-PatternMatched' version exists, the script will also load that data
% e.g. fname_slice_data='Mg 1 Specimen 1 SLICE_0156 SLICE_0156 161';
fname_slice_data='Mg 1 Specimen 1 SLICE_0156 SLICE_0156 161';

%location of MTEX 5.8.1
%if MTEX loaded this doesn't matter
%if MTEX is not loaded, then it will start MTEX up
mtex_location='C:\Users\benja\OneDrive\Documents\MATLAB\mtex-5.8.1'; %working with 5.8.1

[ebsd_original,dataset_header,ebsd_patternmatched] = load_h5oina_pm(fname_slice_data,fname_folder,mtex_location);


%% Plot data
setMTEXpref('xAxisDirection','east');           %aztec
setMTEXpref('zAxisDirection','outofPlane');     %aztec

figure; plot(ebsd_original,ebsd_original.orientations); 
nextAxis
plot(ebsd_patternmatched,ebsd_patternmatched.orientations);