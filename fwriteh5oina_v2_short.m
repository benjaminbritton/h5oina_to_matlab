function fwriteh5oina_v2_short(fileLoc,ebsd,header)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Setup
% Settings.DataName=num2str(Settings.SliceNo);
% Settings.DataName=num2str(1);

%Setup EBSD and Header
dtype=['/EBSD/Data/']; %EBSD data location
htype=['/EBSD/Header/']; %Header data location (e.g. microscope settings)

%change all the map arrays to make dream like them - Type 1 
array_change1=@(x) int32(((x)));

%change all the map arrays to make dream like them - Type 2 
array_change2=@(x) single(((x)));

%change all the map arrays to make dream like them - Type 3
array_change3=@(x) uint8(((x)));

%% General

if exist(fileLoc.HDF5FullFile,'file') == 2 %file exists
    file1s=strfind(fileLoc.HDF5FullFile,'\');
    justfile=fileLoc.HDF5FullFile(file1s(end)+1:end);
    
    text1=input(['The file ' justfile ' exists, overwrite [y/N]?'],'s');
    if strcmpi(text1,'y')
        feval('delete',fileLoc.HDF5FullFile); %#ok<FVAL> 
    else
        error('The file exists, and the h5 writer will not work - choose a different file output name')
    end
end

    if exist(fileLoc.HDF5FullFile,'file') == 0 % do this bit only once
        %NOTE: not reading values here so could manually update or adapt loadEBSD_h5oina_1.m
    
        %Format Version
        data_name=['/' 'Format Version'];
        C='5.0';
        hdf5write(fileLoc.HDF5FullFile,data_name,C);

        %Index
        data_name=['/' 'Index'];
        c=array_change1(1);
        hdf5write(fileLoc.HDF5FullFile,data_name,c,'writemode','append');
        %sort attributes
        h5writeatt( fileLoc.HDF5FullFile ,data_name , 'Type' , 'Single' );

    else
    
    end

%% Data [needs editing]

    %Band Contrast
    m=array_change3(ebsd.prop.Band_Contrast);
    data_name=['/' fileLoc.DataName dtype 'Band Contrast'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');    
    
    %Band Slope
    m=array_change3(ebsd.prop.Band_Slope);
    data_name=['/' fileLoc.DataName dtype 'Band Slope'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');

    %Bands
    m=array_change3(ebsd.prop.Bands);
    data_name=['/' fileLoc.DataName dtype 'Bands'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');

    %Beam Position X [A]
    m=array_change2(ebsd.prop.Beam_Position_X);
    data_name=['/' fileLoc.DataName dtype 'Beam Position X'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');
    h5writeatt( fileLoc.HDF5FullFile ,data_name , 'Unit' , 'um' );

    %Beam Position Y [A]
    m=array_change2(ebsd.prop.Beam_Position_Y);
    data_name=['/' fileLoc.DataName dtype 'Beam Position Y'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');
    h5writeatt( fileLoc.HDF5FullFile ,data_name , 'Unit' , 'um' );

    %Detector Distance
    m=array_change2(ebsd.prop.Detector_Distance);
    data_name=['/' fileLoc.DataName dtype 'Detector Distance'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');

    %Euler [A] - 
%     ebsd.prop.eulerCombo=[ebsd.prop.euler1, ebsd.prop.euler2, ebsd.prop.euler3];
    ebsd.prop.eulerCombo=[ebsd.prop.euler1;ebsd.prop.euler2;ebsd.prop.euler3];
    m=array_change2(ebsd.prop.eulerCombo);
%     m=array_change2(ebsd.orientations.Euler);
%     ebsd.prop.eulerCombo=[ebsd.orientations.phi1; ebsd.orientations.Phi; ebsd.orientations.phi2];
%     m=array_change2(ebsd.prop.eulerCombo);
    m=round(m,4); %ask for help 
    data_name=['/' fileLoc.DataName dtype 'Euler'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');
    h5writeatt( fileLoc.HDF5FullFile ,data_name , 'Unit' , 'rad' );

    %Mean Angular Deviation [A]
    m=array_change2(ebsd.prop.Mean_Angular_Deviation);
    data_name=['/' fileLoc.DataName dtype 'Mean Angular Deviation'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append'); 
    h5writeatt( fileLoc.HDF5FullFile ,data_name , 'Unit' , 'rad' );

    %Pattern Center X
    m=array_change2(ebsd.prop.Pattern_Center_X);
    data_name=['/' fileLoc.DataName dtype 'Pattern Center X'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');

    %Pattern Center Y
    m=array_change2(ebsd.prop.Pattern_Center_Y);
    data_name=['/' fileLoc.DataName dtype 'Pattern Center Y'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');

    %Pattern Quality
    m=array_change2(ebsd.prop.Pattern_Quality);
    data_name=['/' fileLoc.DataName dtype 'Pattern Quality'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');

    %Phase
    m=array_change3(ebsd.prop.Phase);
%         m=array_change3(ebsd.phase); %use this normally
    data_name=['/' fileLoc.DataName dtype 'Phase'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');

    %X [A]
    if isfield(ebsd.prop,'OI_X') == 1
        m=array_change2(ebsd.prop.OI_X); %capital - and playing nice with the updated loader
    else
        m=array_change2(ebsd.prop.X); %capital - and playing nice with the updated loader
    end
%     m=array_change2(ebsd.prop.x);
    data_name=['/' fileLoc.DataName dtype 'X'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');
        h5writeatt( fileLoc.HDF5FullFile ,data_name , 'Unit' , 'um' );

    %Y [A]
      if isfield(ebsd.prop,'OI_Y') == 1
        m=array_change2(ebsd.prop.OI_Y); %capital - and playing nice with the updated loader
    else
        m=array_change2(ebsd.prop.Y); %capital - and playing nice with the updated loader
      end
%     m=array_change2(ebsd.prop.y);
    data_name=['/' fileLoc.DataName dtype 'Y'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');
    h5writeatt( fileLoc.HDF5FullFile ,data_name , 'Unit' , 'um' );


%% Header
  


    %Phases
    htype1=[htype 'Phases/']; 
    %Note: set up for Steel - could automate again
    
        %Color m
        C=array_change2([ebsd.CS.color(1); ebsd.CS.color(2);ebsd.CS.color(3)]');
        data_name=['/' fileLoc.DataName htype1 num2str(1) '/Color'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');

        %Lattice Angles [A] m
        C=array_change2([ebsd.CS.alpha; ebsd.CS.beta; ebsd.CS.gamma]);
        data_name=['/' fileLoc.DataName htype1 num2str(1) '/Lattice Angles'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
            %sort attributes
            h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'rad' );

        %Lattice Dimensions [A] m
        C=array_change2([norm(ebsd.CS.axes(1)); norm(ebsd.CS.axes(2));norm(ebsd.CS.axes(3))]);
        data_name=['/' fileLoc.DataName htype1 num2str(1) '/Lattice Dimensions'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
            %sort attributes
            h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'angstrom' );

        %Laue Group [A] m
        C=array_change2(9); %don't know where this is stored
        data_name=['/' fileLoc.DataName htype1 num2str(1) '/Laue Group'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
            %sort attributes
            h5writeatt( fileLoc.HDF5FullFile , data_name , 'Symbol' , ebsd.CS.LaueName);

        %Number Reflectors m
        C=array_change2(199); %don't know where to find
        data_name=['/' fileLoc.DataName htype1 num2str(1) '/Number Reflectors'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');

        %Phase Id m
        C=array_change2(9);
        data_name=['/' fileLoc.DataName htype1 num2str(1) '/Phase Id'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');

        %Phase Name
        C=ebsd.CS.mineral;
        data_name=['/' fileLoc.DataName htype1 num2str(1) '/Phase Name'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
        
        %Reference
        C='Proc. R. Soc. London, Ser.A[PRLAAZ], vol.A174, page 457';
        data_name=['/' fileLoc.DataName htype1 num2str(1) '/Reference'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
        
        %Space Group [A] m
        C=array_change2(194); %don't know where to find
        data_name=['/' fileLoc.DataName htype1 num2str(1) '/Space Group'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
            %sort attributes
            h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , ebsd.CS.LaueName );              
            
    %Project File
    m=(header.Project_File);
    data_name=['/' fileLoc.DataName htype 'Project File'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');

    %Project Label
    m=(header.Project_Label);
    data_name=['/' fileLoc.DataName htype 'Project Label'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append');
    
    %Scanning Rotation Angle [A] m
    m=array_change2(header.Scanning_Rotation_Angle);
    data_name=['/' fileLoc.DataName htype 'Scanning Rotation Angle'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append'); 
        %sort attributes
        h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'rad' );
    
    %Specimen Orientation Euler [A] m
    m=array_change2(header.Specimen_Orientation_Euler);
    data_name=['/' fileLoc.DataName htype 'Specimen Orientation Euler'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append'); 
        %sort attributes
        h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'rad' );

    %Stage Position %Add if needed
    htype2=[htype 'Stage Position/']; 
        % Rotation [A] m
        C=array_change2(0.19154935); %don't know where to find
        data_name=['/' fileLoc.DataName htype2 'Rotation'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
            %sort attributes
            h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'rad' );

        % Tilt [A] m
        C=array_change2(0.9599136); %don't know where to find
        data_name=['/' fileLoc.DataName htype2 'Tilt'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
            %sort attributes
            h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'rad' );

        % X [A] m
        C=array_change2(1.76028); %don't know where to find
        data_name=['/' fileLoc.DataName htype2 'X'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
            %sort attributes
            h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'mm' );

        % Y [A] m
        C=array_change2(-63.64423); %don't know where to find
        data_name=['/' fileLoc.DataName htype2 'Y'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
            %sort attributes
            h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'mm' );

        % Z [A] m
        C=array_change2(86.44423); %don't know where to find
        data_name=['/' fileLoc.DataName htype2 'Z'];
        hdf5write(fileLoc.HDF5FullFile,data_name,C,'writemode','append');
            %sort attributes
            h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'mm' );

    %X Cells [A] m
    m=array_change1(header.X_Cells);
    data_name=['/' fileLoc.DataName htype 'X Cells'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append'); 
        %sort attributes
        h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'px' );

    %X Step [A] m
    m=array_change2(header.X_Step);
    data_name=['/' fileLoc.DataName htype 'X Step'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append'); 
        %sort attributes
        h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'um' );
    
    %Y Cells [A] m
    m=array_change1(header.Y_Cells);
    data_name=['/' fileLoc.DataName htype 'Y Cells'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append'); 
        %sort attributes
        h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'px' );
    
    %Y Step [A] m
    m=array_change2(header.Y_Step);
    data_name=['/' fileLoc.DataName htype 'Y Step'];
    hdf5write(fileLoc.HDF5FullFile,data_name,m,'writemode','append'); 
        %sort attributes
        h5writeatt( fileLoc.HDF5FullFile , data_name , 'Unit' , 'um' );

end
