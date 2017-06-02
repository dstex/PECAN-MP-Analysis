%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotInterArrvDist_spirals.m
% Created 06/03/2016 by Dan Stechman
% 
% This script is used to plot distributions of interarrival times
% on a spiral-by-spiral basis. These can then be used to visually 
% determine where the local minimum lies between the peaks of the
% bimodal distribution (assuming the case adheres to such a distribution).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clearvars;

flight = '20150706';
probe = 'PIP';

saveFigs = 1;
noDisp = 0;

savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');

proc2File = ['/Users/danstechman/GoogleDrive/PECAN-Data/mp-data/' flight '/proc2.' flight '.' probe '.cdf'];

time_secs = nc_varget(proc2File,'Time_in_seconds');
inter_arr_all=[0;diff(time_secs)];

bins = logspace(-7, 0, 70);

if saveFigs
    saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end
	if (exist([saveDir '/' probe '-IntArrv-Freq-Time'], 'dir') ~= 7)
		mkdir([saveDir '/' probe '-IntArrv-Freq-Time'])
	end
end

for ix=1:length(startT)
    spiralLocs = find(time_secs >= startT(ix) & time_secs < endT(ix));
    intArr = inter_arr_all(spiralLocs);
    
    if saveFigs && noDisp
        figure('visible','off','Position', [10,10,1500,1000]);
    else
        figure('Position', [10,10,1500,1000]);
    end
    histogram(intArr,bins);
    title([flight ' - Spiral ' num2str(ix) ' - ' probe ' Int Arrv Time Freq']);
    set(gca,'Xscale','log');
    ylabel('Frequency');
    xlabel('Interarrival Time [sec]');
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    grid
    if saveFigs
        print([saveDir '/' probe '-IntArrv-Freq-Time/' flight '_' probe '_IntArrvFreq_S' num2str(ix)],'-dpng','-r0')
    end
end