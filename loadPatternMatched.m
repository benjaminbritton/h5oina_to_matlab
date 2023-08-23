%% script portion to load pattern matched data
% load h5oina file after pattern matching
% Note: this is WIP - functional but will be rewritten soon
% 23/08/23 RB
%% MTEX

location_mtex='H:\MTEX\mtex-5.8.1';             % MTEX folder location
run(fullfile(location_mtex,'startup_mtex.m'));  % load MTEX

setMTEXpref('xAxisDirection','east');           %aztec
setMTEXpref('zAxisDirection','outofPlane');     %aztec

%% File locations
fullfile1 = 'C:\Users\rmb07\OneDrive - UBC\7 - Projects\7 - Mg slices\WBV - shared\Mg 1 Specimen 1 SLICE_0156 SLICE_0156 161.h5oina'; %location of original h5oina file
fullfile2 = 'C:\Users\rmb07\OneDrive - UBC\7 - Projects\7 - Mg slices\WBV - shared\Mg 1 Specimen 1 SLICE_0156 SLICE_0156 161-PatternMatching.h5oina'; %location of pattern matched h5oina file

%% Load data - pattern matched

% Load data - pattern matched
[ebsd_o,~,opt_o] = loadEBSD_h5oina_getopt(fullfile1);
[ebsd1,ebsd,~,~,MSdata2] = loadEBSD_h5oina_PatMatch(fullfile2,opt_o) ;

% sort extra parts
ebsd.prop.CrossCorr=MSdata2.Cross_Correlation_Coefficient;
ebsd.prop.Band_Contrast=ebsd1.prop.Band_Contrast;
ebsd.prop.Band_Slope=ebsd1.prop.Band_Slope;
ebsd.prop.Pattern_Quality=ebsd_o.prop.Pattern_Quality;
ebsd.prop.X=ebsd1.prop.X;
ebsd.prop.Y=ebsd1.prop.Y;

% create header
header=ebsd_o.opt.Header;                                   

% optional clean up 
clear ebsd_o opt_o ebsd1 MSdata2  
clear fullfile1 fullfile2 location_mtex

%% Load data - h5oina (normal output)
% simple version if required

%file location
% fullfile1 = 'C:\Users\rmb07\OneDrive - UBC\7 - Projects\7 - Mg slices\WBV - shared\Mg 1 Specimen 1 SLICE_0156 SLICE_0156 161.h5oina'; %location of original h5oina file

% Load data - simple
% [ebsd,~]=loadEBSD_h5oina_1(fullfile1);    % load ebsd
% header=ebsd.opt.Header;                   % get header data

% optional clean up 
% clear fullfile1  location_mtex