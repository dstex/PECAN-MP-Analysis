% imageview - Plot a 2D image record onto the given axis
%     axis_handle - handle pointing to axis you wish to draw the 2D image
%     on
%     imgcdf - netCDF file handle for uncompressed image netCDF
%     frame - # of record to draw
%     rec_num - parent_rec_num from autoanalysis netCDF
%     reject - image_auto_reject from autoanalysis netCDF
function dump_Images(filename, probe)
	FlagColorMap = [1 1 1; 0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 1 1];
	%colormap(FlagColorMap);
	fh = figure('units','inches','outerposition',[0, 0, 8.5, 11],'papertype','usletter');
	pos = get(fh,'Position');
	set(fh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
	
	colormap(FlagColorMap);
	
	rec_num = 888;
	
	%rec_num = 3;
	reject = 0;
	imgcdf = fopen(filename,'r','l');
	
	for frame = 1:20000
		% for frame = 8700:8701
		
		%seekstatus=fseek(imgcdf,frame*4112,'bof');
		%if seekstatus ~= 0
		%    ferror(imgcdf)
		%end
		[date, data]=readChuck(imgcdf);
		hhmmss = date(4)*10000+date(5)*100+date(6);
		millisec = date(7);
		imagedata=ImgDecode(data);
		dims = size(imagedata);
		newimage = zeros(dims(1), dims(2)*8);
		timeentry = date(4)*10000 + date(5)*100 + date(6);
		
		i = 1;
		k = 2;
		y = 1;
		endslice = [170 170 170 170 170 170 170 170];
		endslice1 = [85 0 0 0];
		invalidslice = [-1 -1 -1 -1 -1 -1 -1 -1];
		
		NormalColor = 2;                                               % Black
		FlagStuckBitColor = 3;                                         % Blue
		FlagRejectColor = 3;                                           % Red
		FlagHollowColor = 4;                                           % Green
		
		%particlepos = find(rec_num == frame,1);
		
		CurColor = NormalColor;
		
		% Find the first autoanalysis entry associated with this image
		p = 0;
		
		
		while(k<dims(1))
			
			switch rec_num
				case 222
					if(isequal(imagedata(k,:),endslice1) == 0)
						if(isequal(imagedata(k, :),invalidslice)==0)
							for(i=1:dims(2))
								for(j=1:8)
									if(imagedata(k, i) == -1)
										newimage(y, (i-1)*8+j) = FlagRejectColor;
									else
										newimage(y, (i-1)*8+j) = CurColor-CurColor*bitget(uint8(imagedata(k, i)), 9-j);
									end
								end
							end
							y = y + 1;
							p = p + 1;
						else
							k = dims(1);
						end
						k = k + 1;
					else
						if(y>1)
							newimage(y-1,:)=5;
						end
						newimage(y,:) = 6;
						y = y + 1;
						k = k + 1;
						p = 0;
					end
					
				otherwise
					if(isequal(imagedata(k,:),endslice) == 0 )
						if(isequal(imagedata(k, :),invalidslice)==0)
							
							for(i=1:dims(2))
								for(j=1:8)
									if(imagedata(k, i) == -1)
										newimage(y, (i-1)*8+j) = FlagRejectColor;
									else
										newimage(y, (i-1)*8+j) = CurColor-CurColor*bitget(uint8(imagedata(k, i)), 9-j);
									end
								end
							end
							y = y + 1;
							p = p + 1;
						else
							k = dims(1);
						end
						k = k + 1;
					else
						
						newimage(y,:) = 6;
						y = y + 1;
						k = k + 2;
						p = 0;
						%             if(~isempty(particlepos))
						%                 particlepos = particlepos + 1;
						%                 CurColor = NormalColor;
						%             end
					end
					
			end
			% Get particle header first
			%particleno = uint8(imagedata(k, 1:2));
			%date = uint8(imagedata(k, 3:7));
			%slicecount = bitshift(uint8(imagedata(k,8)), -1);
			%slicecount = mod(uint8(imagedata(k,8)), 128)
			
			
			
			
		end
		
		if y<800
			y=800;
		end
		
		ni = newimage(1:y,:);
		
		% Set aspect ratio to scale - change width to adjust to height of ImgAxes
		
		% title(axes_handle, [num2str(netcdf.getVar(imgcdf,netcdf.inqVarID(imgcdf,'hour'),frame-1, 1)) ':' ...
		%     num2str(netcdf.getVar(imgcdf,netcdf.inqVarID(imgcdf,'minute'),frame-1, 1)) ':' ...
		%     num2str(netcdf.getVar(imgcdf,netcdf.inqVarID(imgcdf,'second'),frame-1, 1)) '.' ...
		%     num2str(netcdf.getVar(imgcdf,netcdf.inqVarID(imgcdf,'millisec'),frame-1, 1))])
		%title(axes_handle, [num2str(imgcdf{'hour'}(frame)) ':' num2str(imgcdf{'minute'}(frame)) ':' num2str(imgcdf{'second'}(frame)) '.' num2str(imgcdf{'millisec'}(frame))])
		iSubplot=mod(frame,20);
		if iSubplot==1
			%    figure(fh)
			%    colormap(FlagColorMap);
			figure(fh)
			pos = get(gca, 'Position');
			pos(1) = 0.03;
			pos(3) = 0.97;
			set(gca, 'Position', pos)
		end
		
		if iSubplot ==0
			iSubplot = 20;
		end
		%colormap(FlagColorMap);
		subplot(20,1,iSubplot)
		imagesc(ni', [1 6]);
		ylabel(['\fontsize{4}' num2str(hhmmss) '.' num2str(millisec)],'FontWeight','normal')
		%axis off;
		axis equal
		box on;
		set(gca,'XTickLabel','')
		set(gca,'YTickLabel','')
		set(gca,'xtick',[])
		set(gca,'ytick',[])
		
		if mod(frame,20)==0
			%frm = getframe( fh );
			print(fh,sprintf( [probe '_%06d.pdf'], hhmmss ),'-dpdf','-r0');
			%imwrite( frm.cdata, sprintf( 'CIP_%06d.png', hhmmss ) ,'XResolution',1800,'YResolution',1800);
			clf
		end
		
		
		
	end
end


function dResult=ImgDecode(data)
	
	bytes=dec2hex(data,2);
	
	i=1;
	ii=1;
	b1full=dec2bin(hex2dec(bytes(:,:)),8);
	b2 = bin2dec(b1full(:,4:8));
	
	while i<4096 & ii<13600
		b1 = b1full(i,:);
		curi = i;
		i=i+1;
		if b1(3) == '1'
			%             i=i+1;
		elseif b1(1) == '0' & b1(2) == '0'
			%            b2=bin2dec(b1(4:8));
			for k=1:b2(curi)+1;
				if i < length(bytes)
					decomp(ii,:)=bytes(i,:);
				else break
				end
				ii=ii+1;
				i=i+1;
			end
		elseif b1(1) == '1' & b1(2) == '0'
			%            b2=bin2dec(b1(4:8));
			for k=1:b2(curi)+1;
				decomp(ii,:)='00';
				ii=ii+1;
			end
		elseif b1(2) == '1' & b1(1) == '0'
			%            b2=bin2dec(b1(4:8));
			for k=1:b2(curi)+1;
				decomp(ii,:)='FF';
				ii=ii+1;
			end
		else
			%kk;
		end
	end
	
	found = 0;
	i=1;
	count=0;
	while found == 0
		if decomp(i,:)=='AA'
			count=count+1;
		else
			count=0;
		end
		
		if count == 8
			found=1;
			dd=i+1:8:length(decomp)-7;
			nWierd=0;
		end
		
		if i==length(decomp)-8 % Add to avoid no 'AA' even though wierd to have no 'AA'...
			found =1;
			nWierd=1;
			nWierdTotal =nWierdTotal +1;
		end
		i=i+1;
	end
	
	
	if nWierd ==0
		decomp_convert=[hex2dec(decomp(dd+7,:)),hex2dec(decomp(dd+6,:)),hex2dec(decomp(dd+5,:)),hex2dec(decomp(dd+4,:)),...
			hex2dec(decomp(dd+3,:)),hex2dec(decomp(dd+2,:)),hex2dec(decomp(dd+1,:)),hex2dec(decomp(dd,:))];
		k2=[decomp(dd,:),decomp(dd+1,:),decomp(dd+2,:),decomp(dd+3,:),decomp(dd+4,:),decomp(dd+5,:),decomp(dd+6,:),decomp(dd+7,:)];
		
		if length(decomp_convert) < 1700
			decomp_convert(length(decomp_convert):1700,:)=-1;
		end
	end
	
	dResult=decomp_convert;
	
end

function [date,data]=readChuck(imgf)
	
	date=fread(imgf,8,'uint16');
	data=fread(imgf,4096,'uint8');
	
end
