clear all;
close all;

config.codepath = 'D:\107_KinSoftChallenge\PostFRET_ExpData\postFRETcodes\';
config.datapath = 'D:\107_KinSoftChallenge\PostFRET_ExpData\expSet3\';
config.resultpath = 'D:\107_KinSoftChallenge\PostFRET_ExpData\expSet3_Results\';



%% Window to select file%%%%

cd(config.datapath);
[config.filelist, path] = uigetfile('*.dat','Choose FRET Files: trace','MultiSelect','on');
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

config.FRETUL = 1.1; % FRET upper limit
config.FRETLL = -0.1; % FRET lower limit


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
    config.datatimestep = data(2,1)- data(1,1);
    d(:,1) = (1:bleachtime)'*config.datatimestep;
    close all;
    fclose 'all'; 
        
    for filenum = 2:size(filenames,2)
        d = load(filenames{filenum});
        datatimestep = d(2,1)-d(1,1);
        
        if datatimestep - config.datatimestep > 1e-10
            display('time step is different for this data set, modify the code to change the time');
            break                    
        end
        d(:,1) = data(end,1)+  (1:size(d,1))'*config.datatimestep; % redo the time:  
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


config.bleachtime = bleachtime*config.datatimestep;

toc

% data format:
% clm  1     2             3            4         5
%  %t (s)	Idd (a.u.)	Ida (a.u.)	Iaa (a.u.)	FRET E  

% datalength:
% accumulative length of each molecule. molecule 1 start with row 1 and
% stops at row datalength(1); molecule 2 starts at datalength(1)+1 and ends
% at row datalength(2)...
data(:,4) = data(:,2) + data(:,3);
a = data(:,3)./data(:,4);
a(isnan(a)) = 0;
a(a>config.FRETUL) = 1; a(a<config.FRETLL) = 0;
data(:,5) = a;

rawdata = data;
cd(config.resultpath);
save('rawdata.mat', 'rawdata', 'datalength', 'config');
cd(config.codepath);

figure; plot(rawdata(:,2)); hold on;
plot(rawdata(:,3)); plot(rawdata(:,4)); title('raw data A D and total counts');
figure; plot(rawdata(:,5)); title('raw data FRET');
