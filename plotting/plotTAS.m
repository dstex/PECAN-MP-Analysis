clearvars;

flight = '20150706';

savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data';

saveFigs	= 1;
noDisp		= 0;
Ftype		= '-dpdf';
% Ftype		= '-dpng';

pipTASf = [dataPath '/mp-data/' flight '/00PIP20150706003806.csv'];
cipTASf = [dataPath '/mp-data/' flight '/01CIP20150706003806.csv'];
casTASf = [dataPath '/mp-data/' flight '/04CAS20150706003806.csv'];
% cdpTASf = [dataPath '/mp-data/' flight '/02CDP20150706003806.csv'];

flRawFile = nc_attget([dataPath '/' flight '_PECANparams.nc'],nc_global,'FL_rawFile');

fltLvlFile = [dataPath '/FlightLevelData/' flRawFile];

loadPIPcsv
loadCIPcsv
loadCAScsv
% loadCDPcsv

FL_True_Air_Speed = nc_varget(fltLvlFile,'TAS.d');
FL_ADDU_True_Air_Speed = nc_varget(fltLvlFile,'TasADDU.1');
FL_HH = nc_varget(fltLvlFile,'HH');
FL_MM = nc_varget(fltLvlFile,'MM');
FL_SS = nc_varget(fltLvlFile,'SS');
FL_Time = (FL_HH*3600) + (FL_MM*60) + FL_SS;
clearvars FL_HH FL_MM FL_SS

PIP_Time_flr = floor(PIP_Time);
CIP_Time_flr = floor(CIP_Time);
CAS_Time_flr = floor(CAS_Time);

gap = diff(CAS_Time_flr);
gapIx = find(gap > 1);
gapSize = gap(gapIx);

if ~isempty(gapIx)
	CAS_newTime = [CAS_Time_flr(1:gapIx(1));(CAS_Time_flr(gapIx(1))+1:CAS_Time_flr(gapIx(1))+gapSize(1)-1)'];
	CAS_newTAS = [CAS_True_Air_Speed(1:gapIx(1));NaN(length(CAS_Time_flr(gapIx(1))+1:CAS_Time_flr(gapIx(1))+gapSize(1)-1),1)];
	
	if length(gapIx) > 2
		for ix = 2:length(gapIx)
			CAS_newTime = [CAS_newTime; CAS_Time_flr(gapIx(ix-1)+1:gapIx(ix)); (CAS_Time_flr(gapIx(ix))+1:CAS_Time_flr(gapIx(ix))+gapSize(ix)-1)'];
			CAS_newTAS = [CAS_newTAS; CAS_True_Air_Speed(gapIx(ix-1)+1:gapIx(ix)); NaN(length(CAS_Time_flr(gapIx(ix))+1:CAS_Time_flr(gapIx(ix))+gapSize(ix)-1),1)];
		end
	end
	
	CAS_newTime = [CAS_newTime; CAS_Time_flr(gapIx(ix)+1:end)];
	CAS_newTAS = [CAS_newTAS; CAS_True_Air_Speed(gapIx(ix)+1:end)];
end

% Determine which for which time range we have data from all probes (for direct comparison)
time1 = max([min(PIP_Time_flr) min(CIP_Time_flr) min(CAS_newTime)]);
time2 = min([max(PIP_Time_flr) max(CIP_Time_flr) max(CAS_newTime)]);

PIPstrt = find(PIP_Time_flr == time1);
PIPend = find(PIP_Time_flr == time2);
CIPstrt = find(CIP_Time_flr == time1);
CIPend = find(CIP_Time_flr == time2);
CASstrt = find(CAS_newTime == time1);
CASend = find(CAS_newTime == time2);
FLstrt = find(FL_Time == time1);
FLend = find(FL_Time == time2);

pipTAS = PIP_True_Air_Speed(PIPstrt:PIPend);
cipTAS = CIP_True_Air_Speed(CIPstrt:CIPend);
casTAS = CAS_newTAS(CASstrt:CASend);
flTAS = FL_True_Air_Speed(FLstrt:FLend);
flADDUTAS = FL_True_Air_Speed(FLstrt:FLend);

if saveFigs && noDisp
	figure('visible','off','Position', [10,10,900,900]);
else
	figure('Position', [10,10,900,900]);
end
scatter(pipTAS,cipTAS,'c');
title('20150706 - PIP TAS vs. CIP TAS','FontSize',18);
xlabel('PIP TAS','FontSize',16);
ylabel('CIP TAS','FontSize',16);
grid
axis equal
refline(1,0)
if saveFigs
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([savePath flight '/FL-Scatters/' flight '_PIP-CIP-TAS_Scatter'],Ftype,'-r0')
end

if saveFigs && noDisp
	figure('visible','off','Position', [10,10,900,900]);
else
	figure('Position', [10,10,900,900]);
end
scatter(pipTAS,casTAS,'c');
title('20150706 - PIP TAS vs. CAS TAS','FontSize',18);
xlabel('PIP TAS','FontSize',16);
ylabel('CAS TAS','FontSize',16);
grid
axis equal
refline(1,0)
if saveFigs
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([savePath flight '/FL-Scatters/' flight '_PIP-CAS-TAS_Scatter'],Ftype,'-r0')
end

if saveFigs && noDisp
	figure('visible','off','Position', [10,10,900,900]);
else
	figure('Position', [10,10,900,900]);
end
scatter(pipTAS,flTAS,'c');
title('20150706 - PIP TAS vs. FL TAS','FontSize',18);
xlabel('PIP TAS','FontSize',16);
ylabel('FL TAS','FontSize',16);
grid
axis equal
refline(1,0)
if saveFigs
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([savePath flight '/FL-Scatters/' flight '_PIP-FL-TAS_Scatter'],Ftype,'-r0')
end

if saveFigs && noDisp
	figure('visible','off','Position', [10,10,900,900]);
else
	figure('Position', [10,10,900,900]);
end
scatter(pipTAS,flADDUTAS,'c');
title('20150706 - PIP TAS vs. FL ADDU TAS','FontSize',18);
xlabel('PIP TAS','FontSize',16);
ylabel('FL ADDU TAS','FontSize',16);
grid
axis equal
refline(1,0)
if saveFigs
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([savePath flight '/FL-Scatters/' flight '_PIP-FLADDU-TAS_Scatter'],Ftype,'-r0')
end
%%

if saveFigs && noDisp
	figure('visible','off','Position', [10,10,900,900]);
else
	figure('Position', [10,10,900,900]);
end
scatter(cipTAS,casTAS,'c');
title('20150706 - CIP TAS vs. CAS TAS','FontSize',18);
xlabel('CIP TAS','FontSize',16);
ylabel('CAS TAS','FontSize',16);
grid
axis equal
refline(1,0)
if saveFigs
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([savePath flight '/FL-Scatters/' flight '_CIP-CAS-TAS_Scatter'],Ftype,'-r0')
end

if saveFigs && noDisp
	figure('visible','off','Position', [10,10,900,900]);
else
	figure('Position', [10,10,900,900]);
end
scatter(cipTAS,flTAS,'c');
title('20150706 - CIP TAS vs. FL TAS','FontSize',18);
xlabel('CIP TAS','FontSize',16);
ylabel('FL TAS','FontSize',16);
grid
axis equal
refline(1,0)
if saveFigs
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([savePath flight '/FL-Scatters/' flight '_CIP-FL-TAS_Scatter'],Ftype,'-r0')
end

if saveFigs && noDisp
	figure('visible','off','Position', [10,10,900,900]);
else
	figure('Position', [10,10,900,900]);
end
scatter(cipTAS,flADDUTAS,'c');
title('20150706 - CIP TAS vs. FL ADDU TAS','FontSize',18);
xlabel('CIP TAS','FontSize',16);
ylabel('FL ADDU TAS','FontSize',16);
grid
axis equal
refline(1,0)
if saveFigs
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([savePath flight '/FL-Scatters/' flight '_CIP-FLADDU-TAS_Scatter'],Ftype,'-r0')
end

%%

if saveFigs && noDisp
	figure('visible','off','Position', [10,10,900,900]);
else
	figure('Position', [10,10,900,900]);
end
scatter(casTAS,flTAS,'c');
title('20150706 - CAS TAS vs. FL TAS','FontSize',18);
xlabel('CAS TAS','FontSize',16);
ylabel('FL TAS','FontSize',16);
grid
axis equal
refline(1,0)
if saveFigs
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([savePath flight '/FL-Scatters/' flight '_CAS-FL-TAS_Scatter'],Ftype,'-r0')
end

if saveFigs && noDisp
	figure('visible','off','Position', [10,10,900,900]);
else
	figure('Position', [10,10,900,900]);
end
scatter(casTAS,flADDUTAS,'c');
title('20150706 - CAS TAS vs. FL ADDU TAS','FontSize',18);
xlabel('CAS TAS','FontSize',16);
ylabel('FL ADDU TAS','FontSize',16);
grid
axis equal
refline(1,0)
if saveFigs
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([savePath flight '/FL-Scatters/' flight '_CAS-FLADDU-TAS_Scatter'],Ftype,'-r0')
end

%%

if saveFigs && noDisp
	figure('visible','off','Position', [10,10,900,900]);
else
	figure('Position', [10,10,900,900]);
end
scatter(flTAS,flADDUTAS,'c');
title('20150706 - FL TAS vs. FL ADDU TAS','FontSize',18);
xlabel('FL TAS','FontSize',16);
ylabel('FL ADDU TAS','FontSize',16);
grid
axis equal
refline(1,0)
if saveFigs
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([savePath flight '/FL-Scatters/' flight '_FL-FLADDU-TAS_Scatter'],Ftype,'-r0')
end