function rgb = HomeyerRainbow(nc)
	path1 = linspace(0.8*pi,1.8*pi,nc);
	path2 = linspace(-0.33*pi,0.33*pi,nc);
	
	y = [linspace(0.3,0.85,nc*2/5) linspace(0.9,0.0,nc - nc*2/5)];
	
	u = 0.4*sin(path1);
	v = 0.55*sin(path2) + 0.1;
	
	rgb_from_yuv = [1, 0, 1.13983; 1, -0.39465, -0.58060; 1, 2.03211, 0];
	
	rgb = NaN(length(y),3);
	
	for i = 1:length(y)
		yuv = [y(i), u(i), v(i)];
% 		rgb(i,:) = rgb_from_yuv*yuv';
% 		rgb(i,:) = yuv*rgb_from_yuv;
		rgb(i,:) = colorspace('YUV->RGB',yuv);
	end
	
% 	rgb(rgb<0) = 0;
% 	rgb(rgb>1) = 1;
end