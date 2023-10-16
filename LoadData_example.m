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
astro_location='C:\Users\benja\OneDrive\Documents\Github\AstroEBSD_v2'; %working with 5.8.1

% run(fullfile(InputUser.Astro_loc,"start_AstroEBSD.m"));
RTM.Phase_Folder = fullfile(astro_location,'phases'); %location of the AstroEBSD phases super-folder


[ebsd_original,dataset_header,ebsd_patternmatched,h5_original,h5_patternmatch] = load_h5oina_pm(fname_slice_data,fname_folder,mtex_location,astro_location);


%% Plot data
setMTEXpref('xAxisDirection','east');           %aztec
setMTEXpref('zAxisDirection','outofPlane');     %aztec

figure; plot(ebsd_original,ebsd_original.orientations); 
nextAxis
plot(ebsd_patternmatched,ebsd_patternmatched.orientations);


figure; plot(ebsd_original,ebsd_original.prop.X); 
nextAxis
plot(ebsd_patternmatched,ebsd_patternmatched.prop.X);

%% Load a pattern
% figure; plot(ebsd_patternmatched,ebsd_patternmatched.prop.CrossCorr); 

%select a point from the map
figure;
plot(ebsd_patternmatched,ebsd_patternmatched.prop.CrossCorr);
colorbar;
[xi,yi]=ginput(1);

%find the point closest to this
r_xy=(ebsd_patternmatched.prop.x-xi).^2+(ebsd_patternmatched.prop.y-yi).^2;
[minr,pattern_number]=min(r_xy);
hold on;
scatter(ebsd_patternmatched(pattern_number).prop.x,ebsd_patternmatched(pattern_number).prop.y);

%%
%extract the pattern
[Pat_Exp]=loadEBSP_T(h5_original,dataset_header,pattern_number);
figure;
imagesc(Pat_Exp);
colormap('gray');
axis equal; axis tight; axis xy;

% imwrite(flipud(uint8(EBSDpat)),'PatTest.png');


%% Sort out the pattern centre convention etc.
PC_pattern_OI=[ebsd_patternmatched(pattern_number).prop.Pattern_Center_X, ebsd_patternmatched(pattern_number).prop.Pattern_Center_Y ebsd_patternmatched(pattern_number).prop.Detector_Distance];
[PC_pattern_Astro,PatternInfo] = PC_OI_to_Astro(PC_pattern_OI,dataset_header);
% PC_pattern_Astro=[0.5 0.5 0.5]
[EBSD_geom ] = EBSP_Gnom( PatternInfo,PC_pattern_Astro); %you can change PC_in if you want

figure;
subplot(1,1,1);
I1=pPattern(Pat_Exp,EBSD_geom);

%% Simulate the pattern


InputUser.Phase_Input  = {'Mg'};

RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},RTM.Phase_Folder);
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

%% simulate using Bruker orientation
% eangs=[221,58,114]*pi/180; %read from a pattern match in dynamics
eangs=[308,82,45]*pi/180;
Detector_tilt=eye(3);
gmatrix=RTM.Rz(eangs(3))*RTM.Rx(eangs(2))*RTM.Rz(eangs(1));
[ Pat_sim_eang ] = EBSP_gen( EBSD_geom,gmatrix*Detector_tilt,screen_int); %generate the EBSP for this iteration 

figure;
subplot(2,1,1);
I1=pPattern(Pat_Exp,EBSD_geom);
subplot(2,1,2);
I1=pPattern(Pat_sim_eang,EBSD_geom);


%%
Settings_Cor=struct;

gmatrix_mtex=ebsd_patternmatched(pattern_number).orientations.matrix;
% gmatrix_mtex=eye(3);
detector_Euler=double(dataset_header.Detector_Orientation_Euler);
g_camera=RTM.Rz(detector_Euler(3))*RTM.Rx(detector_Euler(2))*RTM.Rz(detector_Euler(1));
g_tilt=RTM.Rx(double(dataset_header.Tilt_Angle));

[Pat_Sim]=EBSP_gen( EBSD_geom,gmatrix_mtex'*g_tilt*g_camera',screen_int);

[ Pat_Sim_cor,Settings_Cor_sim] = EBSP_BGCor( Pat_Sim,Settings_Cor );
% Clean up the experimental pattern
Settings_Cor_exp=Settings_Cor;
Settings_Cor_exp.MaskOut=1; %deal with Mask points
Settings_Cor_exp.maskthresh=2; %threshold value for points in the EBSP below this to be set to noise
[ Pat_Exp_cor,Settings_Cor_exp] = EBSP_BGCor( Pat_Exp,Settings_Cor_exp );

figure;
subplot(1,3,1);
pPattern(Pat_Sim_cor,EBSD_geom)
subplot(1,3,2);
pPattern(Pat_Exp_cor,EBSD_geom)

subplot(1,3,3);
pPattern(Pat_Exp_cor-Pat_Sim_cor,EBSD_geom);

%% Exhaustive check of options
%{
G_2=g_tilt*g_camera;
G_3=g_tilt'*g_camera;
G_4=g_tilt'*g_camera';
G_5=g_tilt*g_camera';
G_6=g_camera*g_tilt;
G_7=g_camera'*g_tilt;
G_8=g_camera'*g_tilt';
G_9=g_camera*g_tilt';


[Pat_2]=EBSP_gen( EBSD_geom,gmatrix_mtex'*G_2,screen_int);
[Pat_3]=EBSP_gen( EBSD_geom,gmatrix_mtex'*G_3,screen_int);
[Pat_4]=EBSP_gen( EBSD_geom,gmatrix_mtex'*G_4,screen_int);
[Pat_5]=EBSP_gen( EBSD_geom,gmatrix_mtex'*G_5,screen_int);
[Pat_6]=EBSP_gen( EBSD_geom,gmatrix_mtex'*G_6,screen_int);
[Pat_7]=EBSP_gen( EBSD_geom,gmatrix_mtex'*G_7,screen_int);
[Pat_8]=EBSP_gen( EBSD_geom,gmatrix_mtex'*G_8,screen_int);
[Pat_9]=EBSP_gen( EBSD_geom,gmatrix_mtex'*G_9,screen_int);

figure; 
subplot(3,3,1); pPattern(Pat_Exp,EBSD_geom);
subplot(3,3,2); pPattern(Pat_2,EBSD_geom);title('2');
subplot(3,3,3); pPattern(Pat_3,EBSD_geom);title('3');
subplot(3,3,4); pPattern(Pat_4,EBSD_geom);title('4');
subplot(3,3,5); pPattern(Pat_5,EBSD_geom);title('5');
subplot(3,3,6); pPattern(Pat_6,EBSD_geom);title('6');
subplot(3,3,7); pPattern(Pat_7,EBSD_geom);title('7');
subplot(3,3,8); pPattern(Pat_8,EBSD_geom);title('8');
subplot(3,3,9); pPattern(Pat_9,EBSD_geom);title('9');

%%
figure; 
subplot(3,3,1); pPattern(Pat_Exp,EBSD_geom);
subplot(3,3,2); pPattern(Pat_5,EBSD_geom);title('5');
subplot(3,3,3); pPattern(Pat_7,EBSD_geom);title('7');
subplot(3,3,4); pPattern(Pat_5,EBSD_geom);title('5');
% subplot(3,3,5); pPattern(Pat_5,EBSD_geom);title('5');
% subplot(3,3,6); pPattern(Pat_6,EBSD_geom);title('6');
subplot(3,3,7); pPattern(Pat_7,EBSD_geom);title('7');
% subplot(3,3,8); pPattern(Pat_8,EBSD_geom);title('8');
% subplot(3,3,9); pPattern(Pat_9,EBSD_geom);title('9');
% figure;
% 
%}

%% Pattern match this pattern and see what we get

%setttings for RTM
RTM.screensize = 88; %size of the library patterns and the resize of the raw EBSPs
RTM.Sampling_Freq=8; %Set the SO(3) sampling freq in degrees
RTM.iterations = 4;%Set the number of iterations to do in the refinement step
RTM.LPTsize = 72; %LPT size used in pixels

%From AstroEBSD
%background correction
Settings_CorX.gfilt=0; %use a high pass filter (do you mean high pass?)
Settings_CorX.gfilt_s=2; %low pass filter sigma
Settings_CorX.radius=0; %use a radius mask
Settings_CorX.radius_frac=0.85; %fraction of the pattern width to use as the mask
Settings_CorX.hotpixel=0; %hot pixel correction
Settings_CorX.hot_thresh=1000; %hot pixel threshold
Settings_CorX.resize=1; %resize correction
Settings_CorX.size=RTM.screensize; %image height
Settings_CorX.RealBG=0; %use a real BG
Settings_CorX.EBSP_bgnum=30; %number of real pattern to use for BG
Settings_CorX.SquareCrop = 1; %make square the EBSP
Settings_CorX.SplitBG=0; %deal with a split screen

Settings_CorX.MaskOut=1; %deal with Mask points
Settings_CorX.maskthresh=2; %threshold value for points in the EBSP below this to be set to noise

if Settings_CorX.SquareCrop == 1
    PC_square_x= (PC_pattern_Astro(1)*PatternInfo.ScreenWidth-(PatternInfo.ScreenWidth-PatternInfo.ScreenHeight)/2)/(PatternInfo.ScreenHeight);
    PC_square=PC_pattern_Astro;
    PC_square(1)=PC_square_x;
end

[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );

tic
%correct the experiment
[ Pat_exp_cor ] = EBSP_BGCor( Pat_Exp,Settings_CorX);

%prepare the experimental pattern for refinement
[Pat_Ref_r,XCF_data_fill] = refine_prep(Pat_exp_cor,SettingsXCF,RTM);
[EBSD_geom ] = EBSP_Gnom( RTM,PC_square); %you can change PC_in if you want

%refine
G_start=gmatrix_mtex'*g_tilt*g_camera';

[G_Refined,regout_R] = refine5(Pat_Ref_r,EBSD_geom,EBSD_geom.PC,G_start,SettingsXCF,screen_int,RTM_info.isHex,RTM);
toc
%
[Pat_Sim_In]=EBSP_gen( EBSD_geom,G_start,screen_int);
[Pat_Sim_Ref]=EBSP_gen( EBSD_geom,G_Refined,screen_int);

figure;
subplot(1,3,1); pPattern(Pat_exp_cor,EBSD_geom); title('exp');
subplot(1,3,2); pPattern(Pat_Sim_In,EBSD_geom); title('sim in');
subplot(1,3,3); pPattern(Pat_Sim_Ref,EBSD_geom); title('sim ref');


%% Run more than 1 pattern

tic
PC_x_array=ebsd_patternmatched.prop.Pattern_Center_X;
PC_y_array=ebsd_patternmatched.prop.Pattern_Center_Y;
PC_z_array=ebsd_patternmatched.prop.Detector_Distance;

%an odd way becasue the EBSD container doesn't quite play nice
num_pats=size(ebsd_patternmatched);
num_pats=num_pats(1)*num_pats(2);

refined_G_out=zeros(3,3,num_pats);
refined_R_out=zeros(1,num_pats);

parfor pattern_number=1:num_pats

[Pat_Exp]=loadEBSP_T(h5_original,dataset_header,pattern_number);
%correct the experiment
[ Pat_exp_cor ] = EBSP_BGCor( Pat_Exp,Settings_CorX);

%sort out the PC
PC_pattern_OI=[PC_x_array(pattern_number), PC_y_array(pattern_number) ,PC_z_array(pattern_number)];
[PC_pattern_Astro,PatternInfo] = PC_OI_to_Astro(PC_pattern_OI,dataset_header);

PC_square=PC_pattern_Astro;

if Settings_CorX.SquareCrop == 1
    PC_square_x= (PC_pattern_Astro(1)*PatternInfo.ScreenWidth-(PatternInfo.ScreenWidth-PatternInfo.ScreenHeight)/2)/(PatternInfo.ScreenHeight);
    PC_square(1)=PC_square_x;
end

%prepare the experimental pattern for refinement
[Pat_Ref_r,XCF_data_fill] = refine_prep(Pat_exp_cor,SettingsXCF,RTM);
[EBSD_geom ] = EBSP_Gnom( RTM,PC_square); %you can change PC_in if you want
[G_Refined,regout_R] = refine5(Pat_Ref_r,EBSD_geom,EBSD_geom.PC,G_start,SettingsXCF,screen_int,RTM_info.isHex,RTM);

refined_R_out(:,pattern_number)=regout_R(4);
refined_G_out(:,:,pattern_number)=G_Refined;

%{
%check the answer

[Pat_Sim_In]=EBSP_gen( EBSD_geom,G_start,screen_int);
[Pat_Sim_Ref]=EBSP_gen( EBSD_geom,G_Refined,screen_int);

figure;
subplot(1,3,1); pPattern(Pat_exp_cor,EBSD_geom); title('exp');
subplot(1,3,2); pPattern(Pat_Sim_In,EBSD_geom); title('sim in');
subplot(1,3,3); pPattern(Pat_Sim_Ref,EBSD_geom); title('sim ref');
%}
end
toc
save data1
%%


