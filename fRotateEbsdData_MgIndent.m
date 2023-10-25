function [ebsd]=fRotateEbsdData_MgIndent(ebsd)



%% Calc - can be adjusted

%rotation1 = 53.3degrees about z (sample frame)
% ebsdrot1=rotate(ebsd,53.3*degree);
ebsdrot1=rotate(ebsd,rotation.byAxisAngle(zvector,53.3*degree));%,'keepXY');
%rotation2= 90degrees about x (sample frame)
ebsdrot2=rotate(ebsdrot1,rotation.byAxisAngle(xvector,90*degree),'keepXY'); 
% note: final sample rotation carried out in dream3D pipeline, hence 'keepXY'used

ebsd=ebsdrot2; %output

%% Check plot
% oM=ipfHSVKey(ebsd.CS);
% oM.inversePoleFigureDirection=zvector;
% figure;
% plot(ebsdrot2,oM.orientation2color(ebsdrot2.orientations))