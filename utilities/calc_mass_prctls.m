function [D_10p,D_25p,D_50p,D_75p,D_90p] = calc_mass_prctls(binMid,massSD)
[nTime,nBins] = size(massSD);

D_10p = NaN(nTime,1);
D_25p = NaN(nTime,1);
D_50p = NaN(nTime,1);
D_75p = NaN(nTime,1);
D_90p = NaN(nTime,1);

for t=1:nTime
    
    
    totalMass = nansum(massSD(t,:));
    
    if totalMass == 0 % If total mass for current time is 0 skip iteration
        bRun = 0 ; %Skip to next time
    end
    
    binNum = 1;
    bRun = 1;
    while bRun == 1 && binNum < nBins
        if ( nansum(massSD(t,1:binNum)) <= totalMass*0.10 ) && ( nansum(massSD(t,1:binNum+1)) >= totalMass*0.10 )
            D_10p(t)=(binMid(binNum)+binMid(binNum+1))/2;
            bRun =0 ;
        end
        binNum = binNum+1;
    end
    
    binNum = 1;
    bRun = 1;
    while bRun == 1 && binNum < nBins
        if ( nansum(massSD(t,1:binNum)) <= totalMass*0.25 ) && ( nansum(massSD(t,1:binNum+1)) >= totalMass*0.25 )
            D_25p(t)=(binMid(binNum)+binMid(binNum+1))/2;
            bRun =0 ;
        end
        binNum = binNum+1;
    end
    
    binNum = 1;
    bRun = 1;
    while bRun == 1 && binNum < nBins
        if ( nansum(massSD(t,1:binNum)) <= totalMass*0.5 ) && ( nansum(massSD(t,1:binNum+1)) >= totalMass*0.5 )
            D_50p(t)=(binMid(binNum)+binMid(binNum+1))/2;
            bRun =0 ;
        end
        binNum = binNum+1;
    end
    
    binNum = 1;
    bRun = 1;
    while bRun == 1 && binNum < nBins
        if ( nansum(massSD(t,1:binNum)) <= totalMass*0.75 ) && ( nansum(massSD(t,1:binNum+1)) >= totalMass*0.75 )
            D_75p(t)=(binMid(binNum)+binMid(binNum+1))/2;
            bRun =0 ;
        end
        binNum = binNum+1;
    end
    
    binNum = 1;
    bRun = 1;
    while bRun == 1 && binNum < nBins
        if ( nansum(massSD(t,1:binNum)) <= totalMass*0.90 ) && ( nansum(massSD(t,1:binNum+1)) >= totalMass*0.90 )
            D_90p(t)=(binMid(binNum)+binMid(binNum+1))/2;
            bRun =0 ;
        end
        binNum = binNum+1;
    end
end