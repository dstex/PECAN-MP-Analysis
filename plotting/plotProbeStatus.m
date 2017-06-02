function plotProbeStatus
	clear all; close all;


	savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
	dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

	flight = '20150706';
	
	saveFigs = 1;
	noDisp = 0;

	plotDiodeVoltsCIP	= 1;
	plotLaserCurrCIP	= 1;
	plotDiodeVoltsPIP	= 1;
	plotLaserCurrPIP	= 1;
	
	addTemp = 1;

	if saveFigs
		saveDir = [savePath flight];
		if (exist(saveDir, 'dir') ~= 7)
			mkdir(saveDir)
		end
		if (exist([saveDir '/ProbeStatusPlots'], 'dir') ~= 7 )
			mkdir([saveDir '/ProbeStatusPlots'])
		end
	end
	
	startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
	endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');
	PIP_acptStartT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_acptStartT');
	PIP_acptEndT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_acptEndT');
	PIP_rjctStartT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_rjctStartT');
	PIP_rjctEndT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_rjctEndT');
	
	startT = startT/3600/24;
	endT = endT/3600/24;
	PIP_acptStartT = PIP_acptStartT/3600/24;
	PIP_acptEndT = PIP_acptEndT/3600/24;
	PIP_rjctStartT = PIP_rjctStartT/3600/24;
	PIP_rjctEndT = PIP_rjctEndT/3600/24;
	
	if addTemp
		fltLvlFile = ['/Users/danstechman/GoogleDrive/PECAN-Data/FlightLevelData/Processed/' flight '_FltLvl_Processed.mat'];
		flData = load(fltLvlFile);
		TA = flData.TA;
		timeSecsFL = flData.time_secs_FL;
		timeFLserial = timeSecsFL/3600/24;
	end

	%% Case Selection
	switch flight
		case '20150617'
			csvCIPStr = '01CIP20150617005743.csv';
			csvPIPStr = '00PIP20150617005743.csv';
		case '20150620'
			csvCIPStr = '01CIP20150620005009.csv';
			csvPIPStr = '00PIP20150620005009.csv';
		case '20150701'
			csvCIPStr = '01CIP20150701034344.csv';
			csvPIPStr = '00PIP20150701034344.csv';
		case '20150702'
			csvCIPStr = '01CIP20150702015858.csv';
			csvPIPStr = '00PIP20150702015858.csv';
		case '20150706'
			csvCIPStr = '01CIP20150706003806.csv';
			csvPIPStr = '00PIP20150706003806.csv';
		case '20150709'        
			csvCIPStr = '01CIP20150709002356.csv';
			csvPIPStr = '00PIP20150709002356.csv';
	end

	csvCIPFile = ['/Users/danstechman/GoogleDrive/PECAN-Data/mp-data/' flight '/' csvCIPStr];
	csvPIPFile = ['/Users/danstechman/GoogleDrive/PECAN-Data/mp-data/' flight '/' csvPIPStr];


	%% Extract data from CIP csv status file

	tasfilename = csvCIPFile;
	loadTASinfo
	TimeC = Time;
	timehhmmssC = insec2hhmmss(TimeC);
	Diode_1_VoltsC = Diode_1_Volts;
	Diode_32_VoltsC = Diode_32_Volts;
	Diode_64_VoltsC = Diode_64_Volts;
	Laser_CurrentC = Laser_Current;
	
	serialTimeC = TimeC/3600/24;
	

	%% Create array of different colors for spiral locators
	colors = varycolor(length(startT));
	colors = colors(randperm(length(startT)),:);
	
	
	%% CIP Plot creation

	if plotDiodeVoltsCIP
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,800]);
		else
			figure('Position', [10,10,1200,800]);
		end
		if addTemp
			yyaxis left
		end
		h1 = plot(serialTimeC,Diode_1_VoltsC,'b');
		hold on
		h2 = plot(serialTimeC,Diode_32_VoltsC,'r');
		h3 = plot(serialTimeC,Diode_64_VoltsC,'k');
		h4 = plot(xlim,[1.5 1.5],'b--');
		plot(xlim,[3.2 3.2],'b--');
		h5 = plot(xlim,[2 2],'r--');
		plot(xlim,[3.6 3.6],'r--');
		datetickzoom('x','HH:MM:SS');
		set(gca,'XTickLabelRotation',45);
		title([flight ' - CIP Diode Voltages']);
		xlabel('Time (UTC)');
		ylabel('Diode Voltage (volts)');
		grid
		if addTemp
			yyaxis right
			plot(timeFLserial,TA,'g')
			ylabel(sprintf('Temperature (%cC)',char(176)));
			set(gca,'YDir','reverse');
		end

		for ix=1:length(startT)
			if mod(ix,2) ~= 0
				line([startT(ix) startT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3);
				line([endT(ix) endT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3);
			else
				line([startT(ix) startT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3,'Linestyle','--');
				line([endT(ix) endT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3,'Linestyle','--');
			end
			
		end

		
		set(findall(gcf,'-property','FontSize'),'FontSize',28)

		legend([h1 h2 h3 h4 h5],{'Diode\_1','Diode\_32','Diode\_64','D1 & D64 Limits',...
			'D32 Limits'},'Fontsize',18);

		dcm_obj = datacursormode(gcf);
		set(dcm_obj,'UpdateFcn',{@myupdatefcn});
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if addTemp
				savefig([saveDir '/ProbeStatusPlots/' flight '_CIP_diodeVolts_temp.fig']);
				print([saveDir '/ProbeStatusPlots/' flight '_CIP_diodeVolts_temp'],'-dpdf','-r0')
			else
				savefig([saveDir '/ProbeStatusPlots/' flight '_CIP_diodeVolts.fig']);
				print([saveDir '/ProbeStatusPlots/' flight '_CIP_diodeVolts'],'-dpdf','-r0')
			end
		end

	end

	if plotLaserCurrCIP
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,800]);
		else
			figure('Position', [10,10,1200,800]);
		end
		
		if addTemp
			yyaxis left
		end
		h1 = plot(serialTimeC,Laser_CurrentC,'b');
		hold on
		h2 = plot(xlim,[60 60],'r--');
		plot(xlim,[120 120],'r--');
		datetickzoom('x','HH:MM:SS');
		set(gca,'XTickLabelRotation',45);
		title([flight ' - CIP Laser Current']);
		xlabel('Time (UTC)');
		ylabel('Laser Current (mA)');
		grid
		if addTemp
			yyaxis right
			plot(timeFLserial,TA,'g')
			ylabel(sprintf('Temperature (%cC)',char(176)));
			set(gca,'YDir','reverse');
		end
		
		for ix=1:length(startT)
			if mod(ix,2) ~= 0
				line([startT(ix) startT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3);
				line([endT(ix) endT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3);
			else
				line([startT(ix) startT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3,'Linestyle','--');
				line([endT(ix) endT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3,'Linestyle','--');
			end
		end

		
		set(findall(gcf,'-property','FontSize'),'FontSize',28)

		legend([h1 h2],{'Laser Current','Current Limits'},'Fontsize',18);

		dcm_obj = datacursormode(gcf);
		set(dcm_obj,'UpdateFcn',{@myupdatefcn});
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if addTemp
				savefig([saveDir '/ProbeStatusPlots/' flight '_CIP_laserCurrent_temp.fig']);
				print([saveDir '/ProbeStatusPlots/' flight '_CIP_laserCurrent_temp'],'-dpdf','-r0')
			else
				savefig([saveDir '/ProbeStatusPlots/' flight '_CIP_laserCurrent.fig']);
				print([saveDir '/ProbeStatusPlots/' flight '_CIP_laserCurrent'],'-dpdf','-r0')
			end
		end
		
	end

	%% Extract data from PIP csv status file

	tasfilename = csvPIPFile;
	loadTASinfo
	
	TimeP = Time;
	timehhmmssP = insec2hhmmss(TimeP);
	Diode_1_VoltsP = Diode_1_Volts;
	Diode_32_VoltsP = Diode_32_Volts;
	Diode_64_VoltsP = Diode_64_Volts;
	Laser_CurrentP = Laser_Current;
	
	serialTimeP = TimeP/3600/24;


	%% PIP Plot creation

	if plotDiodeVoltsPIP
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,800]);
		else
			figure('Position', [10,10,1200,800]);
		end
		
		if addTemp
			yyaxis left
		end
		h1 = plot(serialTimeP,Diode_1_VoltsP,'b');
		hold on
		h2 = plot(serialTimeP,Diode_32_VoltsP,'r');
		h3 = plot(serialTimeP,Diode_64_VoltsP,'k');
		h4 = plot(xlim,[1.5 1.5],'b--');
		plot(xlim,[3 3],'b--');
		h5 = plot(xlim,[1.2 1.2],'r--');
		plot(xlim,[3.7 3.7],'r--');
		datetickzoom('x','HH:MM:SS');
		set(gca,'XTickLabelRotation',45);
		title([flight ' - PIP Diode Voltages']);
		xlabel('Time (UTC)');
		ylabel('Diode Voltage (volts)');
		grid
		if addTemp
			yyaxis right
			plot(timeFLserial,TA,'g')
			ylabel(sprintf('Temperature (%cC)',char(176)));
			set(gca,'YDir','reverse');
		end
		
		
		for ix=1:length(startT)
			if mod(ix,2) ~= 0
				line([startT(ix) startT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3);
				line([endT(ix) endT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3);
			else
				line([startT(ix) startT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3,'Linestyle','--');
				line([endT(ix) endT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3,'Linestyle','--');
			end
			
			if ~isnan(PIP_rjctStartT(ix))
				line([PIP_rjctStartT(ix) PIP_rjctStartT(ix)],get(gca,'YLim'),'Color','k','Linewidth',2,'Linestyle','-.');
			end
		end

		
		set(findall(gcf,'-property','FontSize'),'FontSize',28)

		legend([h1 h2 h3 h4 h5],{'Diode\_1','Diode\_32','Diode\_64','D1 & D64 Limits',...
			'D32 Limits'},'Fontsize',18);


		dcm_obj = datacursormode(gcf);
		set(dcm_obj,'UpdateFcn',{@myupdatefcn});
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if addTemp
				savefig([saveDir '/ProbeStatusPlots/' flight '_PIP_diodeVolts_temp.fig']);
				print([saveDir '/ProbeStatusPlots/' flight '_PIP_diodeVolts_temp'],'-dpdf','-r0')
			else
				savefig([saveDir '/ProbeStatusPlots/' flight '_PIP_diodeVolts.fig']);
				print([saveDir '/ProbeStatusPlots/' flight '_PIP_diodeVolts'],'-dpdf','-r0')
			end
		end
		
	end

	if plotLaserCurrPIP
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,800]);
		else
			figure('Position', [10,10,1200,800]);
		end
		
		if addTemp
			yyaxis left
		end
		h1 = plot(serialTimeP,Laser_CurrentP,'b');
		hold on
		h2 = plot(xlim,[60 60],'r--');
		plot(xlim,[120 120],'r--');
		datetickzoom('x','HH:MM:SS');
		set(gca,'XTickLabelRotation',45);
		title([flight ' - PIP Laser Current']);
		xlabel('Time (UTC)');
		ylabel('Laser Current (mA)');
		grid
		if addTemp
			yyaxis right
			plot(timeFLserial,TA,'g')
			ylabel(sprintf('Temperature (%cC)',char(176)));
			set(gca,'YDir','reverse');
		end
			
		for ix=1:length(startT)
			if mod(ix,2) ~= 0
				line([startT(ix) startT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3);
				line([endT(ix) endT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3);
			else
				line([startT(ix) startT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3,'Linestyle','--');
				line([endT(ix) endT(ix)],get(gca,'YLim'),'Color',colors(ix,:),'Linewidth',3,'Linestyle','--');
			end
			
			if ~isnan(PIP_rjctStartT(ix))
				line([PIP_rjctStartT(ix) PIP_rjctStartT(ix)],get(gca,'YLim'),'Color','k','Linewidth',2,'Linestyle','-.');
			end
		end
		
		set(findall(gcf,'-property','FontSize'),'FontSize',28)

		legend([h1 h2],{'Laser Current','Current Limits'},'Fontsize',18);

		dcm_obj = datacursormode(gcf);
		set(dcm_obj,'UpdateFcn',{@myupdatefcn});
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if addTemp
				savefig([saveDir '/ProbeStatusPlots/' flight '_PIP_laserCurrent_temp.fig']);
				print([saveDir '/ProbeStatusPlots/' flight '_PIP_laserCurrent_temp'],'-dpdf','-r0')
			else
				savefig([saveDir '/ProbeStatusPlots/' flight '_PIP_laserCurrent.fig']);
				print([saveDir '/ProbeStatusPlots/' flight '_PIP_laserCurrent'],'-dpdf','-r0')
			end
		end
		
	end
end

% Used to modify datatip for PECAN probe status plots
function txt = myupdatefcn(~,event_obj)
    % Customizes text of data tips
    pos = get(event_obj,'Position');
    I = get(event_obj, 'DataIndex');
% 	txt = {['Time: ',num2str(insec2hhmmss(pos(1)*3600*24))],...
% 			['Y: ',num2str(pos(2))]};
	txt = {['Time: ',datestr(pos(1),'HH:MM:SS.FFF')],...
			['Y: ',num2str(pos(2))]};
    alldatacursors = findall(gcf,'type','hggroup');
    set(alldatacursors,'FontSize',20);
end