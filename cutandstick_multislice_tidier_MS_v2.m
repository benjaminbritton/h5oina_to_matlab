%% Cut & Stick New Version


%% Admin - File locations

%location of MTEX 5.8.1
%if MTEX loaded this doesn't matter
%if MTEX is not loaded, then it will start MTEX up
% mtex_location='C:\Users\benja\OneDrive\Documents\MATLAB\mtex-5.8.1'; %working with 5.8.1
mtex_location='H:\MTEX\mtex-5.8.1'; %working with 5.8.1

% Astro location
astro_location='C:\Users\rmb07\Documents\GitHub\AstroEBSD_v2';

% run(fullfile(InputUser.Astro_loc,"start_AstroEBSD.m"));
RTM.Phase_Folder = fullfile(astro_location,'phases'); %location of the AstroEBSD phases super-folder

% Output file name & destination
OutputFolder='E:\Mg_voilume\analysed\MATLAB outputs\';

%% Setup

%start mtex if needed
try EBSD;
catch
    run(fullfile(mtex_location,"startup.m"));
end

%start astro if needed
try astro_loadcheck
catch
    run(fullfile(astro_location,"start_AstroEBSD.m"));
end


%% Choices

% mtex prefs
setMTEXpref('xAxisDirection','east');       %aztec
setMTEXpref('zAxisDirection','outofPlane'); %aztec

% step size/ voxel size (um)
ss=0.15;

% y height for slice [this becomes z height]
y_No=5;

%slice numbers
startSliceNo=81;    %start
finishSliceNo=254;  %finish

% Output file
NameO=['MS_Slices_' num2str(startSliceNo) '-' num2str(finishSliceNo) '_height_' num2str(y_No) '_new.h5oina'];
fout=struct;
fout.HDF5FullFile=[OutputFolder NameO];
fout.DataName='1'; 

%% Loop to load and select data (CUT)

%Get file names list
fname_folder = 'E:\Mg_voilume\analysed\';    %folder
f_List = dir(fname_folder);                        % get list of files
fileNames = {f_List(~[f_List.isdir]).name};         % get file names
fileNames = strrep(fileNames, '.h5oina', '');       % remove file type
getpm = endsWith(fileNames,"-PatternMatching");     % remove duplicates (pattern matched) from list
fileNames=fileNames(~getpm);                        % final list

num_slices=length(fileNames);

for loopNo=1:num_slices

    % file name for slice
    fname_slice_data=fileNames{length(fileNames)+1-loopNo};

    % load data
    [ebsd_original,dataset_header,ebsd_patternmatched,h5_original,h5_patternmatch] = load_h5oina_pm(fname_slice_data,fname_folder,mtex_location,astro_location);

%     figure; plot(ebsd_patternmatched,ebsd_patternmatched.orientations);

    % rotate the data
    [ebsd_patternmatched]=fRotateEbsdData_MgIndent(ebsd_patternmatched);

%     figure; plot(ebsd_patternmatched,ebsd_patternmatched.orientations);  
%     figure; plot(ebsd_patternmatched,ebsd_patternmatched.prop.Band_Contrast);  

    % sort the ebsd -> grid
    ebsd_patternmatched=ebsd_patternmatched.gridify;

    % pick the points for the slice through [TO BE IMPROVED]
    %x1=[min(ebsd.x) max(ebsd.x)];
    x1=[-19.95 30];
    % x1=[-13.95 24]; %lower (y=-16)
    % x1=[-13.95 12.9]; %lower (y=-32)
    y1=[-y_No  -y_No];

    % define the line and cut
    lineSecN=[x1(1) y1(1); x1(2) y1(2)];
    ebsd_lineN=spatialProfile(ebsd_patternmatched,lineSecN);

    %fix a potential issue - mainly for unrotated data
    ebsd_lineN=ebsd_lineN(ebsd_lineN.y==max(ebsd_lineN.y)); %if it picks two lines
    
    % find the number of points at the desired step size
    ebsd_lineN_x=round((x1(2)-x1(1))/ss);

    % make a container for the new strip of ebsd
    ebsd_lineN1=repmat(ebsd_lineN(1),ebsd_lineN_x,1);

    % loop time - points along the line
    for m=1:ebsd_lineN_x

        % x position
        val1=x1(1)+(ss*(m-1)); 
        ind = interp1(ebsd_lineN.x,1:length(ebsd_lineN.x),val1,'nearest','extrap');

        % unindexed points
        if ebsd_lineN(ind).phase==0
            ebsd_lineN1(m).prop.euler1=0.0;       
            ebsd_lineN1(m).prop.euler2=0.0;
            ebsd_lineN1(m).prop.euler3=0.0;
            ebsd_lineN1(m).prop.Phase=0;
            ebsd_lineN1(m).phase=0;
            ebsd_lineN1(m).phaseId=0;

        else
            %indexed points
            ori=ebsd_lineN(ind).orientations;
            ebsd_lineN1(m).prop.euler1=ori.phi1;       
            ebsd_lineN1(m).prop.euler2=ori.Phi;
            ebsd_lineN1(m).prop.euler3=ori.phi2;
            ebsd_lineN1(m).prop.Phase=ebsd_lineN(ind).phase; %to be improved
            ebsd_lineN1(m).phase=ebsd_lineN(ind).phase;
            ebsd_lineN1(m).orientations=ebsd_lineN(ind).orientations;
        end
        %Beam_Position_Y is used in the h5 writer, and y is thrown away.
        %y is used by MTEX to plot (and created from Beam_Position_Y when
        %the h5 is loaded)
        %these two variables must be the same
        ebsd_lineN1(m).prop.Beam_Position_Y=(num_slices-loopNo+1)*ss;
        ebsd_lineN1(m).y=ebsd_lineN1(m).prop.Beam_Position_Y;

        ebsd_lineN1(m).prop.Beam_Position_X=val1;
        ebsd_lineN1(m).x=ebsd_lineN1(m).prop.Beam_Position_X;
%       ebsd_lineN1(m).y=(ss*length(fileNames))-ss*loopNo; %maybe?  
        
        ebsd_lineN1(m).prop.OI_X=ebsd_lineN(ind).prop.OI_X;
        ebsd_lineN1(m).prop.OI_Y=ebsd_lineN(ind).prop.OI_Y; %come back to!
        
        
        ebsd_lineN1(m).prop.id=(loopNo-1)*(ebsd_lineN_x)+m;

        ebsd_lineN1(m).prop.Band_Contrast=ebsd_lineN(ind).prop.Band_Contrast;
        ebsd_lineN1(m).prop.Band_Slope=ebsd_lineN(ind).prop.Band_Slope;
        ebsd_lineN1(m).prop.Bands=ebsd_lineN(ind).prop.Bands;
        ebsd_lineN1(m).prop.Detector_Distance=ebsd_lineN(ind).prop.Detector_Distance;
        ebsd_lineN1(m).prop.Mean_Angular_Deviation=ebsd_lineN(ind).prop.Mean_Angular_Deviation;
        ebsd_lineN1(m).prop.Pattern_Quality=ebsd_lineN(ind).prop.Pattern_Quality;
        ebsd_lineN1(m).prop.Phase=ebsd_lineN(ind).prop.Phase;  
        ebsd_lineN1(m).phaseId=ebsd_lineN(ind).prop.Phase;
        ebsd_lineN1(m).prop.CrossCorr=ebsd_lineN(ind).prop.CrossCorr;
        ebsd_lineN1(m).phase=ebsd_lineN(ind).prop.Phase;

    end
    
    ebsd_lineT{loopNo}=ebsd_lineN1;

end


% put all the ebsd together
ebsd_lineTall=[ebsd_lineT{:}];

% figure;
plot(ebsd_lineTall,ebsd_lineTall.prop.x);

%orientations
% orientations - complicated version
oriList_e1=ebsd_lineT{1}.prop.euler1;
oriList_e2=ebsd_lineT{1}.prop.euler2;
oriList_e3=ebsd_lineT{1}.prop.euler3;

% for n=1:9%length(fileNames1)
for n=2:length(fileNames)
    oriList_e1=[oriList_e1;ebsd_lineT{n}.prop.euler1];
    oriList_e2=[oriList_e2;ebsd_lineT{n}.prop.euler2];
    oriList_e3=[oriList_e3;ebsd_lineT{n}.prop.euler3];
end

cs_Mg=ebsd_patternmatched.CS; %move this

for n=1:length(oriList_e1)
    % turn them back into orientations
    oriList(n)=orientation('euler',oriList_e1(n),oriList_e2(n),oriList_e3(n),cs_Mg,1);
end

oriList(isnan(oriList))=0; % put zeros - this is a fudge but aztec loads the file this way

%% Make new EBSD container (STICK)

%CS list
CSList=ebsd_lineTall.CSList;

%prop
prop.Band_Contrast=ebsd_lineTall.prop.Band_Contrast;
prop.Band_Slope=ebsd_lineTall.prop.Band_Slope;
prop.Bands=ebsd_lineTall.prop.Bands;
prop.Beam_Position_X=ebsd_lineTall.prop.Beam_Position_X;
prop.Beam_Position_Y=ebsd_lineTall.prop.Beam_Position_Y;

prop.x=prop.Beam_Position_X;
prop.y=prop.Beam_Position_Y;

prop.Detector_Distance=ebsd_lineTall.prop.Detector_Distance;
prop.euler1=oriList.phi1(:)';%in radians
prop.euler2=oriList.Phi(:)';%in radians
prop.euler3=oriList.phi2(:)';%in radians
% prop.euler1=oriList_e1; %alt
% prop.euler2=oriList_e2; %alt
% prop.euler3=oriList_e3; %alt
prop.Mean_Angular_Deviation=ebsd_lineTall.prop.Mean_Angular_Deviation;
prop.Pattern_Center_X=ebsd_lineTall.prop.Pattern_Center_X;
prop.Pattern_Center_Y=ebsd_lineTall.prop.Pattern_Center_Y;
prop.Pattern_Quality=ebsd_lineTall.prop.Pattern_Quality;
prop.Phase=ebsd_lineTall.prop.Phase;
prop.OI_X = ebsd_lineTall.prop.OI_X;
prop.OI_Y= ebsd_lineTall.prop.OI_Y;
% prop.x=ebsd_lineTall.prop.x;%ebsd_lineTall.x;
% prop.y=sort(ebsd_lineTall.prop.y,'descend');%ebsd_lineTall.y;
prop.id=ebsd_lineTall.id;
prop.CrossCorr=ebsd_lineTall.prop.CrossCorr;

ebsd_lineTall.prop.Phase(isnan(ebsd_lineTall.prop.Phase))=0;

% create new EBSD container
ebsd_test1=EBSD(oriList,ebsd_lineTall.prop.Phase,CSList,prop,'unitcell',ebsd_lineTall.unitCell);
% ebsd_test1=EBSD(oriList,ebsd_lineTall.prop.Phase,CSList,prop);

% plot to check
figure; plot(ebsd_test1,ebsd_test1.orientations);
figure; plot(ebsd_test1,ebsd_test1.prop.Band_Contrast);

%sort header [mainly using the last loaded one]
header=dataset_header;
header.X_Cells=ebsd_lineN_x; % no of points in x (pixels)
header.Y_Cells=length(fileNames); % no of points in y (pixels)
header.X_Step=ss; % step size
header.Y_Step=ss; % step size

%% Output - write the h5oina file

ebsd_test1_new=ebsd_test1;

%if you need to rewrite the position data to meddle, you can do here:
% ebsd_test1_new.prop.Beam_Position_Y=ss*loopNo-ebsd_test1_new.prop.Beam_Position_Y;
% ebsd_test1_new.prop.y=ebsd_test1_new.prop.Beam_Position_Y;
% ebsd_test1_new.prop.x=ebsd_test1_new.prop.Beam_Position_X;

fwriteh5oina_v2_short(fout,ebsd_test1_new,header);


%% Read the h5oina file back in & plot

[ebsd_re,dataset_header_re,ebsd_patternmatched_re,h5_original_re,h5_patternmatch_re] = load_h5oina_pm(NameO(1:end-7),OutputFolder,mtex_location,astro_location);

%%

% figure; 
% plot(ebsd_test1_new,ebsd_test1_new.prop.Beam_Position_Y);
% % title('Data before writing to h5oina');
% nextAxis;
% plot(ebsd_test1_new,ebsd_test1_new.prop.y);
% % title('Data after writing and loading')
%%
figure; 
plot(ebsd_test1_new,ebsd_test1_new.prop.Band_Contrast);
title('Data before writing to h5oina');
axis on;
nextAxis;
plot(ebsd_re,ebsd_re.prop.Band_Contrast);
axis on;
title('Data after writing and loading');