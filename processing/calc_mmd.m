% Function to calculate median mass diameter
% Written by Wei Wu

% massSD should be unnormalized by bindwidth prior to being used here (i.e., units of g cm-3)
% iwc should have same units as massSD
% units of bin_mids should correspond to that used in massSD and iwc (most likely cm)

function mmd=calc_mmd(bin_mids,massSD,iwc)

sizeAll=size(massSD);
nTime=sizeAll(1);
nBins=sizeAll(2);
mmd=zeros(nTime,1);

for i=1:nTime
    j=1;
    bRun=1;
    
    if iwc(i) == 0
        mmd(i)=0;
        bRun=0;
	end
	
    while bRun==1 && j<nBins
        if ( nansum(massSD(i,1:j),2)<=iwc(i)/2 ) && ( nansum(massSD(i,1:j+1),2)>=iwc(i)/2 )
            mmd(i)=(bin_mids(j)+bin_mids(j+1))/2;
            bRun=0;
		end
        j=j+1;
    end
end