
% Tool for loading an XML file and running it using the Marta system.
% Run, then select a valid .xml file.
%
% KS 2024/1/5

clear;

parent_directory = uigetdir('',"Select working directory");
cd(parent_directory);

[fName,path] = uigetfile('*.xml',"Select XML file");

[trials,info] = ParseExpFile(fullfile(path,fName));

marta(trials,info);