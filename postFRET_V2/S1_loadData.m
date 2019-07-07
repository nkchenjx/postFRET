% coded by Jixin Chen @ Ohio University,   chenj@ohio.edu   2019


clear all;
close all;

% config.codepath = 'D:\PostFRET_level1\'; % the folder of the codes, change it when in a different path
% config.datapath = 'D:\sim_190212_194543_level1\'; % the folder of the raw data, change it when in a different path
% config.resultpath = 'D:\Results_level1'; % Change it to any folder that you want the results to be saved



%% Window to select file%%%%

% cd(config.datapath);
cd('.\exampleData\');
[config.filelist, path] = uigetfile('*.txt','Choose FRET Files: trace','MultiSelect','on'); % select all traces
PrevFolder=cd;
cd(path);
if iscell(config.filelist) % if multiselection is on, e.filelist is a cell
     nloops = length(config.filelist);
     config.filelist=sort(config.filelist);
elseif ischar(config.filelist) % if multiselection is off, e.filelist is char string
     nloops = 1;
else
     warndlg('e.filelist has not been read properly!')
     return
end
% e.filelist = filelist
fname = config.filelist;


%% load files
tic;
filenames = config.filelist;
data = [];
bleachtime = [];
datalength = [];
if nloops > 1
    d = load(filenames{1});
    data = [d];
    datalength(1) = size(d,1);
    bleachtime = [bleachtime; size(d,1)];
    close all;
    fclose 'all'; 
        
    for filenum = 2:size(filenames,2)
        d = load(filenames{filenum});
        d(:,1) = d(:,1) + data(end,1) + d(2,1)-d(1,1); % redo the time:
        bleachtime = [bleachtime; size(d,1)];
        data = [data; d];
        datalength(filenum,1) = size(d,1)+ datalength(filenum-1,1);
        close all;
        fclose 'all';
    end
    
% elseif nloops == 1
%     d = load(filenames{1});
%     data = [d];
%     datalength(1) = size(d,1);
%     close all;
%     fclose 'all'; 
end

config.datatimestep = data(2,1)- data(1,1);
config.bleachtime = bleachtime*config.datatimestep;

toc

% data format:
% clm  1     2             3            4         5
%  %t (s)	Idd (a.u.)	Ida (a.u.)	Iaa (a.u.)	FRET E  

% datalength:
% accumulative length of each molecule. molecule 1 start with row 1 and
% stops at row datalength(1); molecule 2 starts at datalength(1)+1 and ends
% at row datalength(2)...

rawdata = data;
%cd(config.resultpath);
cd('..\exampleResults');
save('rawdata.mat', 'rawdata', 'datalength', 'config');
%cd(config.codepath);
cd('..');
