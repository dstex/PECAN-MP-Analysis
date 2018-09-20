function plot_CIPfits_inspect
clearvars; close all;

%% Specify various plotting/calculation parameters
flight = '20150620';

avgTime = 10;

zoom = 0;

contLevs = 70;
logCont = 1;

rmvMLmass = 1;

fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_12mm'];


% Standard plots, but using hybrid (CIP obs + CIP extended) SDs

plotNDtemp			= 1;
plotMDtemp			= 1;



if zoom
    diamLim = [0.1 2.5];
else
    diamLim = [0.1 6];
end

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';




%% Load in struct of all original SD data and the CIP Fitted Data
sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.' num2str(avgTime) 'secAvg.mat'];
sDistF = load(sDistFile);

load([dataPath 'mp-data/' flight '/sDist/' flight fileIdStr '.mat']);

loopVctr = 1:length(sprlNames);
% loopVctr = [2];

%% Specify flight-specific plotting parameters
tempRangeAll = [-18.5 22];

if zoom
    NDLim = [1e-7 30];
    NDavgLim = [3e-4 5];
    MDLim = [1e-10 5e-5];
    MDavgLim = [1e-9 1e-5];
else
    NDLim = [1e-7 30];
    NDavgLim = [1e-6 5];
    MDLim = [1e-10 5e-5];
    MDavgLim = [1e-9 1e-5];
end
NDLogCLim = [-5 2];
NDLogLim = [-7 2];
MDLogCLim = [-8 -4];
MDLogLim = [-10 -4];




%% Standard plots
if plotNDtemp
    if any(isinf(diamLim)) || isempty(diamLim)
        diamLim = [min(cipExt_binMid_mm) max(cipExt_binMid_mm)];
    end
    for ix = loopVctr
        
        cipConc = cipConc_cm4_hybrid_igf.(sprlNames{ix});
        cipConc(cipConc == 0) = NaN;
        
        tempCsprl = sDistF.tempC_orig.(sprlNames{ix});
        time_fl = sDistF.time_secsFL_orig.(sprlNames{ix});
        
        time_secs = sDistF.time_secs_avg.(sprlNames{ix});
        alt_avg = sDistF.alt_avg.(sprlNames{ix});
        RH_avg = sDistF.RH_avg.(sprlNames{ix});
        
        if rmvMLmass
            tempCavg = sDistF.tempC_avg.(sprlNames{ix});
            if ~isnan(mlTopTemp(ix))
                Tgt0IXs = find(tempCavg >= mlTopTemp(ix));
            else
                Tgt0IXs = find(tempCavg >= 0);
            end
        end
        
        [D10,D25,D50,D75,D90] = calc_mass_prctls(cipExt_binEdges_mm,cipConc);
        D10(Tgt0IXs) = NaN;
        D25(Tgt0IXs) = NaN;
        D50(Tgt0IXs) = NaN;
        D75(Tgt0IXs) = NaN;
        D90(Tgt0IXs) = NaN;
        
        figure('Position', [10,10,1000,1500]);
        
        
        set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
        
        ax = gca;
        contourf(cipExt_binMid_mm,time_secs/24/3600,log10(cipConc),linspace(NDLogLim(1),NDLogLim(2),contLevs),'LineColor','none');
        hold on;
        %             plot(D10,time_secs/24/3600,'k--','LineWidth',1.5);
        plot(D25,time_secs/24/3600,'k--','LineWidth',3);
        plot(D50,time_secs/24/3600,'k-','LineWidth',3);
        plot(D75,time_secs/24/3600,'k--','LineWidth',3);
        %             plot(D90,time_secs/24/3600,'k--','LineWidth',1.5);
        % 			h = pcolor(cipExt_binMid_mm,time_secs/24/3600,log10(cipConc));
        % 			set(h,'EdgeColor','none')
        xlabel('D (mm)');
        ylabel('Time');
        set(ax,'XMinorTick','on');
        colormap(HomeyerRainbow(contLevs));
        c=colorbar;
        set(c,'Location','southoutside');
        ylabel(c,'log_{10}N(D) (cm^{-4})');
        set(ax, 'CLim', NDLogCLim);
        
        if logCont
            set(ax,'XScale','log');
            labLocFac = 0.83;
            scaleStr = '_log';
            diamLim = [0.15 diamLim(2)];
        else
            labLocFac = 0.95;
            scaleStr = '';
        end
        set(ax,'XLim',diamLim);
        
        % Plot ML top/bottom locations and annotate with temp
        topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
        botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
        hold on
        plot(diamLim, [1 1]*mlTopTime(ix)/24/3600,'k--')
        tMT = text(diamLim(2)*labLocFac,mlTopTime(ix)/24/3600,topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
        plot(diamLim, [1 1]*mlBotTime(ix)/24/3600,'k--')
        tMB = text(diamLim(2)*labLocFac,mlBotTime(ix)/24/3600,botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
        hold off
        
        % Swap time direction depending on spiral direction and define top/bottom spiral temps
        if tempCsprl(1) < tempCsprl(end)
            set(ax,'YDir','reverse');
        end
        
        % Set the tick frequency for the time axis and adjust plot accordingly
        delta = 60; % 1 minute
        dtStrt = (floor(time_secs(1)/delta)*delta)/24/3600;
        dtEnd = time_secs(end)/24/3600;
        dtDelta = delta/24/3600;
        set(ax,'YTick',dtStrt:dtDelta:dtEnd);
        datetick('y','HH:MM','keepticks','keeplimits');
        
        
        title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);
        
        set(findall(gcf,'-property','FontSize'),'FontSize',26)
        set(tMB,'FontSize',14);
        set(tMT,'FontSize',14);
        set(gca,'layer','top'); % Puts tick marks and such on top of the pcolor
        
        dcm_obj = datacursormode(gcf);
        set(dcm_obj,'UpdateFcn',{@myupdatefcn,tempCavg,time_secs,alt_avg,RH_avg});
        
    end
end

if plotMDtemp
    if any(isinf(diamLim)) || isempty(diamLim)
        diamLim = [min(cipExt_binMid_mm) max(cipExt_binMid_mm)];
    end
    for ix = loopVctr
        cipMass = cipMass_gcm4_hybrid_igf.(sprlNames{ix});
        cipMass(cipMass == 0) = NaN;
        time_secs = sDistF.time_secs_avg.(sprlNames{ix});
        
        tempC = sDistF.tempC_orig.(sprlNames{ix});
        time_fl = sDistF.time_secsFL_orig.(sprlNames{ix});
        alt_avg = sDistF.alt_avg.(sprlNames{ix});
        RH_avg = sDistF.RH_avg.(sprlNames{ix});
        
        % Set any mass values in the ML to NaN
        if rmvMLmass
            tempCavg = sDistF.tempC_avg.(sprlNames{ix});
            if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
                mlIXs = find(tempCavg <= mlBotTemp(ix) & tempCavg >= mlTopTemp(ix));
                cipMass(mlIXs,:) = NaN;
                Tgt0IXs = find(tempCavg >= mlTopTemp(ix));
            else
                Tgt0IXs = find(tempCavg >= 0);
            end
        end
        
        [D10,D25,D50,D75,D90] = calc_mass_prctls(cipExt_binEdges_mm,cipMass);
        D10(Tgt0IXs) = NaN;
        D25(Tgt0IXs) = NaN;
        D50(Tgt0IXs) = NaN;
        D75(Tgt0IXs) = NaN;
        D90(Tgt0IXs) = NaN;
        
        
        figure('Position', [10,10,1000,1500]);
        
        
        set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
        
        ax = gca;
        contourf(cipExt_binMid_mm,time_secs/24/3600,log10(cipMass),linspace(MDLogLim(1),MDLogLim(2),contLevs),'LineColor','none');
        hold on;
        plot(D10,time_secs/24/3600,'k--','LineWidth',1.5);
        plot(D25,time_secs/24/3600,'k--','LineWidth',3);
        plot(D50,time_secs/24/3600,'k-','LineWidth',3);
        plot(D75,time_secs/24/3600,'k--','LineWidth',3);
        plot(D90,time_secs/24/3600,'k--','LineWidth',1.5);
        % 			h = pcolor(cipExt_binMid_mm,time_secs/24/3600,log10(mass_twc));
        % 			set(h,'EdgeColor','none')
        xlabel('D (mm)');
        ylabel('Time');
        set(ax,'XMinorTick','on');
        colormap(HomeyerRainbow(contLevs));
        c=colorbar;
        set(c,'Location','southoutside');
        ylabel(c,'log_{10}M(D) (g cm^{-4})');
        set(ax, 'CLim', MDLogCLim);
        if logCont
            set(ax,'XScale','log');
            labLocFac = 0.83;
            scaleStr = '_log';
            diamLim = [0.15 diamLim(2)];
        else
            labLocFac = 0.95;
            scaleStr = '';
        end
        set(ax,'XLim',diamLim);
        
        % Plot ML top/bottom locations and annotate with temp
        topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
        botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
        hold on
        plot(diamLim, [1 1]*mlTopTime(ix)/24/3600,'k--')
        tMT = text(diamLim(2)*labLocFac,mlTopTime(ix)/24/3600,topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
        plot(diamLim, [1 1]*mlBotTime(ix)/24/3600,'k--')
        tMB = text(diamLim(2)*labLocFac,mlBotTime(ix)/24/3600,botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
        hold off
        
        % Swap time direction depending on spiral direction and define top/bottom spiral temps
        if tempC(1) < tempC(end)
            set(ax,'YDir','reverse');
        end
        
        % Set the tick frequency for the time axis and adjust plot accordingly
        delta = 60; % 1 minute
        dtStrt = (floor(time_secs(1)/delta)*delta)/24/3600;
        dtEnd = time_secs(end)/24/3600;
        dtDelta = delta/24/3600;
        set(ax,'YTick',dtStrt:dtDelta:dtEnd);
        datetick('y','HH:MM','keepticks','keeplimits');
        
        
        title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);
        
        set(findall(gcf,'-property','FontSize'),'FontSize',26)
        set(tMB,'FontSize',14);
        set(tMT,'FontSize',14);
        set(gca,'layer','top'); % Puts tick marks and such on top of the pcolor
        
        dcm_obj = datacursormode(gcf);
        set(dcm_obj,'UpdateFcn',{@myupdatefcn,tempCavg,time_secs,alt_avg,RH_avg});
        
    end
end
end

function txt = myupdatefcn(~,event_obj,tempCavg,time_secs,alt_avg,RH_avg)
% Customizes text of data tips
pos = get(event_obj,'Position');
% I = get(event_obj, 'DataIndex');
tempCcall = tempCavg(find(floor(time_secs) == floor(pos(2)*3600*24)));
altcall = alt_avg(find(floor(time_secs) == floor(pos(2)*3600*24)));
rhcall = RH_avg(find(floor(time_secs) == floor(pos(2)*3600*24)));
txt = {sprintf('%.3f mm',pos(1)),...
    datestr(pos(2),'HH:MM:SS'),...
    sprintf('%.3g',10^pos(3)),...
    sprintf('%.2f%cC',tempCcall,char(176)),...
    sprintf('%.2f %%',rhcall),...
    sprintf('%.2f km',altcall/1000)};
alldatacursors = findall(gcf,'type','hggroup');
set(alldatacursors,'FontSize',16);
end