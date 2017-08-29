function [PerOn]=DutyCycleCalculator(Ca); 
sz=length(Ca);
t=1:sz;
% t=st:ed;
%              Ca=medfilt1(Ca);
%             Ca=interp1(t,Ca,tint);
%             fitCa=polyfit(t,Ca,3);
%             corrCa=polyval(fitCa, t);
%             Ca=Ca./corrCa;
  %           plot(Ca)
 
% drawnow
             [pLoc, pAmp]=peakfinder(Ca);
            [tLoc,tAmp]=peakfinder(Ca,[],[],-1);
            amp=mean(pAmp)-mean(tAmp);
            
%            scatter(pLoc,pAmp)
%            scatter(tLoc,tAmp)
%            drawnow

%      

time=0;
count=1;
clear loc val

% figure; hold on;

thresh=1+(1-min(tAmp))/4;

% plot(t,ones(length(t),1)*min(tAmp))
% plot(t, ones(length(t),1)*thresh)
% scatter(t,Ca)

for k=1:length(t)
  
%     amp=(mean(pAmp)+mean(tAmp))/2;
    if Ca(k)>=thresh
%        if Ca(k)>=amp*0.1
%       if Ca(k)>amp
        time=time+1;
        loc(count)=k;
        val(count)=Ca(k);
        count=count+1;
%     end
    end
end

% scatter(loc,val, 'filled')
% drawnow;



PerOn=time/sz;
%scatter(loc,val, 'filled')
%drawnow
%%%%%%%%%%Other option for duty cycle %%%%%%%%%%%%
% for k=1:length(t)
% %     amp=1.0;
%     amp=(mean(pAmp)+mean(tAmp))/2;
% %     if Ca(k)>=amp*0.05+mean(tAmp)
%       if Ca(k)>=amp
%         time=time+1;
%         loc(count)=k;
%         val(count)=Ca(k);
%         count=count+1;
%       end
% end




% scatter(loc,val,'filled')
% drawnow
% pause(2)




   
   