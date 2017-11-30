function FLtempLookup_time(flight,InTimehhmmss)
	% clearvars;
	
	
	dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';
	
	FLfile = [dataPath 'FlightLevelData/Processed/' flight '_FltLvl_Processed.mat'];
	
	load(FLfile, 'TA','time_secs_FL','RH_hybrid');
	
	% InTimehhmmss = input('Time (hhmmss): ');
	
	inTime_secs = hhmmss2insec(InTimehhmmss);
	
	
	[~, inqIxFL] = min(abs(time_secs_FL - inTime_secs));
	
	fprintf('\nNearest FL time (hhmmss): %d\n',insec2hhmmss(time_secs_FL(inqIxFL)));
	fprintf('Nearest FL time (sec): %d\n',time_secs_FL(inqIxFL));
	fprintf('Temp = %f C\n',TA(inqIxFL));
	fprintf('RH = %f%%\n',RH_hybrid(inqIxFL));
	
end