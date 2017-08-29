
figure; hold on;
PerOnMap=zeros(sx,sy);
for i=1:length(x)
    
    dat=Ims(x(i)-5:x(i)+5,y(i)-5:y(i)+5,:);
%     dat=Ims(x(i)-20:x(i)+20,y(i)-20:y(i)+20,:);
    dat=nanmean(dat,2);
    dat=nanmean(dat,1);
    dat=dat(:);
    f=polyfit(t',dat,2);
    feval=polyval(f,t');
    %dat=dat-feval;
    %dat=detrend(dat)';
    dat=dat./polyval(f,t');
    
   if ActMask(x(i),y(i))>0
        PerOnMap(x(i),y(i))=DutyCycleCalculator(dat);
   else
        PerOnMap(x(i), y(i))=0;
   end
   
end

%%
Ims=WT.images;
[sx,sy,sz]=size(Ims);
cell1=createMask(imfreehand);
cell2=createMask(imfreehand);
cell3=createMask(imfreehand);

%%
maskstack1 = repmat(cell1,[1 1 sz]);
maskedStack1 = Ims .* maskstack1;
TC(1,:)=mean(mean(maskedStack1,2),1);

maskstack2 = repmat(cell2,[1 1 sz]);
maskedStack2 = Ims .* maskstack2;
TC(2,:)=mean(mean(maskedStack2,2),1);

maskstack3 = repmat(cell3,[1 1 sz]);
maskedStack3 = Ims .* maskstack3;
TC(3,:)=mean(mean(maskedStack3,2),1);
%%
% sz=170;

% TC(1,:)=WT.CorrCellTCs(1,2:sz);
% TC(2,:)=WT.CorrCellTCs(2,2:sz);
% TC(3,:)=WT.CorrCellTCs(3,2:sz);



t=1:sz;
figure; hold on;
for i=1:3


    dat=TC(i,:);
    f=polyfit(t,dat,2);
    feval=polyval(f,t);
    dat=dat./feval;

    PerOn=DutyCycleCalculator(dat);
end
