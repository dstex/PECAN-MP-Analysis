% This script will create plots comparing size distributions for datasets with
% and without shatter removal and/or shatter reacceptance

function plotSizeDistShatrChecks
close all;clearvars;

flight = '20150706';

probe = 'CIP';

plotND                  = 1;
plotNDtime              = 1;
plotIntArvShatrRccpt    = 1;

saveFigs    = 1;
noDisp      = 1;

savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

sDist_noShatters = ['/Users/danstechman/GoogleDrive/PECAN-Data/mp-data/' flight '/sDist/sdistCI.' flight '.' probe '.noShatters.mat'];
sDist_withShatters = ['/Users/danstechman/GoogleDrive/PECAN-Data/mp-data/' flight '/sDist/sdistCI.' flight '.' probe '.withShatters.mat'];

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');
if strcmp(probe,'CIP')
	intar_threshold = nc_varget([dataPath '/' flight '_PECANparams.nc'],'CIP_intArvThrsh');
else
	intar_threshold = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_intArvThrsh');
end


% Variables from file where shatter removal was applied
importVars = {'timehhmmss','cip2_binmin','cip2_binmid','cip2_conc_minR','cip2_n','bad_cip2_conc_minR','bad_cip2_n',...
    'shatrReject_times','shatrReject_intArr','shatrReject_diam','rccptReject_times','rccptReject_intArr',...
    'rccptReject_diam','loopedTimes','loopedIntArr','loopedDiam'};

load(sDist_noShatters,importVars{:});
time_secs = hhmmss2insec(timehhmmss);
cip2_conc_minR = cip2_conc_minR';
cip2_n = cip2_n';
bad_cip2_conc_minR = bad_cip2_conc_minR';
bad_cip2_n = bad_cip2_n';

% Variables from file where NO shatter removal was applied
importVars_wS = {'cip2_conc_minR','cip2_n'};

temp = load(sDist_withShatters,importVars_wS{:});

cip2_conc_minR_wS = (temp.(importVars_wS{1}))';
cip2_n_wS = temp.(importVars_wS{2});

% cip2_conc_minR_diff = cip2_conc_minR_wS-cip2_conc_minR;



if plotND
    for ix = 1:length(startT)
        timeLocs = find(time_secs >= startT(ix) & time_secs < endT(ix));
        
        
        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,3000,700]);
        else
            figure('Position', [10,10,3000,700]);
        end
        
        subplot(1,3,1)
        stairs(cip2_binmin, nanmean(cip2_conc_minR_wS(:,timeLocs),2), 'r', 'LineWidth', 2);
%         title([flight ' - Spiral ' num2str(ix) ' - CIP WithShatters']);
        title([flight ' - S16 Subset - CIP WithShatters']);
        xlabel('D [mm]');
        ylabel('N(D) [cm^{-4}]');
        set(gca,'Yscale','log');
        ylim([5e-5 130]);
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        grid;
		
        subplot(1,3,2)
        stairs(cip2_binmin, nanmean(cip2_conc_minR(:,timeLocs),2), 'r', 'LineWidth', 2);
%         title([flight ' - Spiral ' num2str(ix) ' - CIP NoShatters']);
        title([flight ' - S16 Subset - CIP NoShatters']);
        xlabel('D [mm]');
        ylabel('N(D) [cm^{-4}]');
        set(gca,'Yscale','log');
        ylim([5e-5 130]);
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        grid;
        
        subplot(1,3,3)
        stairs(cip2_binmin, nanmean(bad_cip2_conc_minR(:,timeLocs),2), 'r', 'LineWidth', 2);
%         title([flight ' - Spiral ' num2str(ix) ' - CIP (Rejected Only)']);
        title([flight ' - S16 Subset - CIP (Rejected Only)']);
        xlabel('D [mm]');
        ylabel('N(D) [cm^{-4}]');
        set(gca,'Yscale','log');
        ylim([5e-5 130]);
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        grid;
        
        tightfig;
        
        if saveFigs
% 			print([savePath flight '/CIP-sizeDist_RejShtrd/' flight '_CIP_sizeDist_rejDfdShtrd_S' num2str(ix)],'-dpng','-r0')
% 			print([savePath flight '/CIP-sizeDist_RejAll/' flight '_CIP_sizeDist_rejDfdAll_S' num2str(ix)],'-dpng','-r0')
			print([savePath flight '/CIP-sizeDist_RejAll/' flight '_CIP_sizeDist_rejDfdAll_S16-Sub'],'-dpng','-r0')
        end
    end
end

if plotNDtime
    for ix = 1:length(startT)
        timeLocs = find(time_secs >= startT(ix) & time_secs < endT(ix));
        
        
        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,3000,700]);
        else
            figure('Position', [10,10,3000,700]);
        end

        subplot(1,3,1)
        contourf(time_secs(timeLocs)/24/3600,cip2_binmid,log10(cip2_conc_minR_wS(:,timeLocs)),-4:0.1:2,'LineColor','none');
%         title([flight ' - Spiral ' num2str(ix) ' - CIP WithShatters']);
        title([flight ' - S16 Subset - CIP WithShatters']);
        ylabel('D [mm]');
        %xlim([startT(ix)/3600/24, endT(ix)/3600/24])
        datetick('x','HH:MM');
        set(gca,'XTickLabelRotation',45);
        colormap(jetmod); %Uses modified 'jet' colormap
        c=colorbar;
        ylabel(c,'log_{10}N(D) [cm^{-4}]');
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        set(gca, 'CLim', [-4,  2]);
        
        subplot(1,3,2)
        contourf(time_secs(timeLocs)/24/3600,cip2_binmid,log10(cip2_conc_minR(:,timeLocs)),-4:0.1:2,'LineColor','none');
%         title([flight ' - Spiral ' num2str(ix) ' - CIP NoShatters']);
        title([flight ' - S16 Subset - CIP NoShatters']);
        ylabel('D [mm]');
        %xlim([startT(ix)/3600/24, endT(ix)/3600/24])
        datetick('x','HH:MM');
        set(gca,'XTickLabelRotation',45);
        colormap(jetmod); %Uses modified 'jet' colormap
        c=colorbar;
        ylabel(c,'log_{10}N(D) [cm^{-4}]');
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        set(gca, 'CLim', [-4,  2]);
        
        subplot(1,3,3)
        contourf(time_secs(timeLocs)/24/3600,cip2_binmid,log10(bad_cip2_conc_minR(:,timeLocs)),-4:0.1:2,'LineColor','none');
%         title([flight ' - Spiral ' num2str(ix) ' - CIP (Rejected Only)']);
        title([flight ' - S16 Subset - CIP (Rejected Only)']);
        ylabel('D [mm]');
        %xlim([startT(ix)/3600/24, endT(ix)/3600/24])
        datetick('x','HH:MM');
        set(gca,'XTickLabelRotation',45);
        colormap(jetmod); %Uses modified 'jet' colormap
        c=colorbar;
        ylabel(c,'log_{10}N(D) [cm^{-4}]');
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        set(gca, 'CLim', [-4,  2]);
        
        tightfig;
        
        if saveFigs
%			print([savePath flight '/CIP-ND-vs-Time_RejShtrd/' flight '_CIP_NDtime_rejDfdShtrd_S' num2str(ix)],'-dpng','-r0')
% 			print([savePath flight '/CIP-ND-vs-Time_RejAll/' flight '_CIP_NDtime_rejDfdAll_S' num2str(ix)],'-dpng','-r0')
			print([savePath flight '/CIP-ND-vs-Time_RejAll/' flight '_CIP_NDtime_rejDfdAll_S16-Sub'],'-dpng','-r0')
        end
    end
end

if plotIntArvShatrRccpt
    for ix = 1:length(startT)
        timeLocs = find(loopedTimes >= startT(ix) & loopedTimes < endT(ix));
        shatrTimeLocs = find(shatrReject_times >= startT(ix) & shatrReject_times < endT(ix));
        rccptTimeLocs = find(rccptReject_times >= startT(ix) & rccptReject_times < endT(ix));
        
        loopedTimes_sub = loopedTimes(timeLocs);
        loopedDiam_sub = loopedDiam(timeLocs);
        
        yPos = intar_threshold(ix);
        
        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,2200,800]);
        else
            figure('Position', [10,10,2200,800]);
        end
        plot(loopedTimes(timeLocs)/3600/24,loopedIntArr(timeLocs),'b.',...
            'MarkerSize',15);
        set(gca,'Yscale','log');
		
		title([flight ' - Spiral ' num2str(ix) ' - CIP Interarrival (shatters/reaccepts Flagged)']);
%         title([flight ' - S16 Subset - CIP Interarrival (shatters/reaccepts Flagged)']);
		
        hold on
        plot(get(gca,'xlim'), [yPos yPos]);
        plot(shatrReject_times(shatrTimeLocs)/3600/24,shatrReject_intArr(shatrTimeLocs),'r.',...
            'MarkerSize',22);
        plot(rccptReject_times(rccptTimeLocs)/3600/24,rccptReject_intArr(rccptTimeLocs),'g.',...
            'MarkerSize',22);
        datetick('x','HH:MM:SS');
        hold off
        ylabel('Interarrival Time [sec]');
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        
        dcm_obj = datacursormode(gcf);
        set(dcm_obj,'UpdateFcn',{@myupdatefcn,loopedTimes_sub,loopedDiam_sub,timeLocs});
        
        if saveFigs
            tightfig;
            print([savePath flight '/CIP-IntArrv_Flags/' flight '_CIP_IntArrv-Flags_S' num2str(ix)],'-dpng','-r0')
% 			print([savePath flight '/CIP-IntArrv_Flags/' flight '_CIP_IntArrv-ShatrRejFlags_S16-Sub'],'-dpng','-r0')
        end
    end
end

end

function txt = myupdatefcn(~,event_obj,loopedTimes,loopedDiam,timeLocs)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
txt = {['Time: ',datestr(loopedTimes(I)/3600/24,'HH:MM:SS.FFF')],...
        ['Ind: ',num2str(timeLocs(I))],...
        ['IntArrv: ',num2str(pos(2))],...
        ['Diam: ',num2str(loopedDiam(I))]};
alldatacursors = findall(gcf,'type','hggroup');
set(alldatacursors,'FontSize',16);
end
