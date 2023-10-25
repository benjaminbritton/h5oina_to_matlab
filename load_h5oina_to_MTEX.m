close all
clear
home
%folder that contains the h5oina files, 
% % e.g. fname_folder='C:\Users\benja\OneDrive\Documents\MATLAB\h5oina';
% fname_folder='C:\Users\benja\OneDrive\Documents\MATLAB\h5oina';
% fname_folder='H:\Ben';
% fname_folder='E:\Mg_voilume\analysed\MATLAB outputs';



%h5oina file, excluding final file extension & folder
%if the '-PatternMatched' version exists, the script will also load that data
% e.g. fname_slice_data='Mg 1 Specimen 1 SLICE_0156 SLICE_0156 161';
% fname_slice_data='MS_Slices_101-120_XZplane_Yheight_16_5';

% fname_folder='E:\SiMaps2\h5oina';
% fname_slice_data='Si_Patterns_32x32 Si StepSize5_160x160um Resolution_SS5';
% phase_plot='Silicon';

% Pick your phase
fname_folder='E:\Mg_voilume\analysed\';
fname_slice_data='Mg 1 Specimen 1 SLICE_0156 SLICE_0156 161';
phase_plot='Magnesium';


%%
% astro_location='C:\Users\benja\OneDrive\Documents\Github\AstroEBSD_v2'; %working with 5.8.1
astro_location='C:\Users\rmb07\Documents\GitHub\AstroEBSD_v2';

%location of MTEX 5.8.1
%if MTEX loaded this doesn't matter
%if MTEX is not loaded, then it will start MTEX up
% mtex_location='C:\Users\benja\OneDrive\Documents\MATLAB\mtex-5.8.1'; %working with 5.8.1
mtex_location='H:\MTEX\mtex-5.8.1'; %working with 5.8.1

%if you are going to pattern match - you need this loadedat
% run(fullfile(InputUser.Astro_loc,"start_AstroEBSD.m"));
RTM.Phase_Folder = fullfile(astro_location,'phases'); %location of the AstroEBSD phases super-folder

[ebsd_original,dataset_header,ebsd_patternmatched,h5_original,h5_patternmatch] = load_h5oina_pm(fname_slice_data,fname_folder,mtex_location,astro_location);

%% Plot data
setMTEXpref('xAxisDirection','east');           %aztec
setMTEXpref('zAxisDirection','outofPlane');     %aztec

figure; plot(ebsd_original,ebsd_original.orientations); 
nextAxis
plot(ebsd_patternmatched,ebsd_patternmatched.orientations);

%%
figure;
plot(ebsd_patternmatched,ebsd_patternmatched.orientations);

oM_H=ipfHSVKey(ebsd_patternmatched(phase_plot).CS);

% IPF
Dirx=xvector;
Diry=yvector;
Dirz=zvector;

figure; 

plot(ebsd_patternmatched(phase_plot),ebsd_patternmatched(phase_plot).prop.Band_Contrast,'micronbar','on');
title('BC');
nextAxis;

oM_H.inversePoleFigureDirection = Dirx; 
plot(ebsd_patternmatched(phase_plot),oM_H.orientation2color(ebsd_patternmatched(phase_plot).orientations),'micronbar','on');
title('X');
nextAxis;

oM_H.inversePoleFigureDirection = Diry; 
plot(ebsd_patternmatched(phase_plot),oM_H.orientation2color(ebsd_patternmatched(phase_plot).orientations),'micronbar','on');
title('Y');
nextAxis;

oM_H.inversePoleFigureDirection = Dirz; 
plot(ebsd_patternmatched(phase_plot),oM_H.orientation2color(ebsd_patternmatched(phase_plot).orientations),'micronbar','on');
title('Z');
nextAxis;

%%
figure;
s1=subplot(3,3,1);
plot(ebsd_patternmatched,ebsd_patternmatched.prop.Beam_Position_X,'micronbar','on','parent',s1);
title('Beam X')
colorbar;

s2=subplot(3,3,2);
plot(ebsd_patternmatched,ebsd_patternmatched.prop.Beam_Position_Y,'micronbar','on','parent',s2);
title('Beam Y')
colorbar;

% s4=subplot(4,3,4);
% plot(ebsd_patternmatched,ebsd_patternmatched.prop.X,'micronbar','on','parent',s4);
% title('X')
% colorbar;
% 
% s5=subplot(4,3,5);
% plot(ebsd_patternmatched,ebsd_patternmatched.prop.Y,'micronbar','on','parent',s5);
% title('Y')
% colorbar;

s7=subplot(3,3,4);
plot(ebsd_patternmatched,ebsd_patternmatched.prop.x,'micronbar','on','parent',s7);
title('position x')
colorbar;

s8=subplot(3,3,5);
plot(ebsd_patternmatched,ebsd_patternmatched.prop.y,'micronbar','on','parent',s8);
title('position y')
colorbar;

s10=subplot(3,3,7);
plot(ebsd_patternmatched,ebsd_patternmatched.prop.Pattern_Center_X,'micronbar','on','parent',s10);
title('PC_X')
colorbar;

s11=subplot(3,3,8);
plot(ebsd_patternmatched,ebsd_patternmatched.prop.Pattern_Center_Y,'micronbar','on','parent',s11);
title('PC_Y')
colorbar;

s12=subplot(3,3,9);
plot(ebsd_patternmatched,ebsd_patternmatched.prop.Detector_Distance,'micronbar','on','parent',s12);
title('DD')
colorbar;

%%

%%
figure;
plot(ebsd_patternmatched,ebsd_patternmatched.prop.Band_Contrast,'micronbar','off');
%%
figure;
s_bc=subplot(3,3,5);
plot(ebsd_patternmatched,ebsd_patternmatched.prop.Band_Contrast,'micronbar','on','parent',s_bc);
title('BC');
hold on;
indexlist(1)=sub2ind([32 32],1,1); %[X_num,Y_num,row,col] - col = x position, row = y position
indexlist(2)=sub2ind([32 32],32,1); %[X_num,Y_num,row,col]
indexlist(3)=sub2ind([32 32],1,32); %[X_num,Y_num,row,col]
indexlist(4)=sub2ind([32 32],32,32); %[X_num,Y_num,row,col]

%make sure the index list is the same size as the colours & locations in
%the plot
ind_loc=[1,3,7,9];
ind_col={'r','g','b','m'};


PC_x_array=ebsd_patternmatched.prop.Pattern_Center_X;
PC_y_array=ebsd_patternmatched.prop.Pattern_Center_Y;
PC_z_array=ebsd_patternmatched.prop.Detector_Distance;

PatData.ScreenWidth=    double(dataset_header.Pattern_Width);
PatData.ScreenHeight=    double(dataset_header.Pattern_Height);

for n=1:numel(indexlist)
    s_pat=subplot(3,3,ind_loc(n));
    pattern_number=indexlist(n);
    PC_pattern_OI=[PC_x_array(pattern_number), PC_y_array(pattern_number) ,PC_z_array(pattern_number)];
    [PC_pattern_Astro,PatternInfo] = PC_OI_to_Astro(PC_pattern_OI,dataset_header);
    [EBSD_geom ] = EBSP_Gnom( PatData,PC_pattern_Astro); %you can change PC_in if you want
    [Pat_Exp]=loadEBSP_T(h5_original,dataset_header,pattern_number);
    pPattern(Pat_Exp,EBSD_geom,s_pat); hold on;
    plot(EBSD_geom.x_screen,EBSD_geom.y_screen(1)*ones(size(EBSD_geom.x_screen)),'Color',ind_col{n},'LineWidth',3);
    plot(EBSD_geom.x_screen,EBSD_geom.y_screen(end)*ones(size(EBSD_geom.x_screen)),'Color',ind_col{n},'LineWidth',3);

    plot(EBSD_geom.x_screen(1)*ones(size(EBSD_geom.y_screen)),EBSD_geom.y_screen,'Color',ind_col{n},'LineWidth',3);
    plot(EBSD_geom.x_screen(end)*ones(size(EBSD_geom.y_screen)),EBSD_geom.y_screen,'Color',ind_col{n},'LineWidth',3);

    mh=mean(EBSD_geom.y_screen);
    plot(EBSD_geom.x_screen,mh*ones(size(EBSD_geom.x_screen)),'Color',ind_col{n},'LineWidth',1,'LineStyle','-');
     mv=mean(EBSD_geom.x_screen);
    plot(mv*ones(size(EBSD_geom.y_screen)),EBSD_geom.y_screen,'Color',ind_col{n},'LineWidth',1,'LineStyle','-');


    scatter(ebsd_patternmatched(indexlist(n)).prop.x,ebsd_patternmatched(indexlist(n)).prop.y,'parent',s_bc,'MarkerEdgeColor',ind_col{n});
end

% indexlist=[1,32,1024-32+1,1024]

%%
odf = calcDensity(ebsd_patternmatched(phase_plot).orientations);
h = [Miller(0,0,0,1,odf.CS),Miller(1,1,-2,0,odf.CS),Miller(1,0,-1,0,odf.CS)];
figure;
plotPDF(odf,h,'antipodal','silent','grid','upper','projection','eangle','grid_res',10*degree)

%% write this to a h5 file
fout_folder=[fname_folder(1:end-1) '_test'];
fout_fname=[fname_slice_data '_out.h5oina'];

fout_full=struct;
fout_full.HDF5FullFile=fullfile(fout_folder,fout_fname);
fout_full.DataName='1'; 
fwriteh5oina_v2_short(fout_full,ebsd_patternmatched,dataset_header);