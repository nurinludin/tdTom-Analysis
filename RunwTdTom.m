function [DataOut]=Run(data, tdTomMask)
%% LOADING IMAGES

if length(fieldnames(data))==1
    mask=[];
    silentcell=[];
    bkgcell=[];
else
    mask=data.mask;
    silentcell=data.silentcell;
    bkgcell=data.bkgcell;
  scthresh=data.SCthresh;
  bkgthresh=data.BKGthresh;
 
end


R=bfopen(data.Location);


pics=R{1};
pics=pics(:,1);


try R{4}.getPlaneDeltaT(0,2)

if ~isempty(R{4}.getPlaneDeltaT(0,0))
for i=1:length(pics)
T(i)=R{4}.getPl/aneDeltaT(0,i-1);
end
else
    T=1:size(pics,3);
end

catch
    
% for i=1:length(pics)
% T(i)=R{4}.getPlaneDeltaT(i-1,0);
% end    
    
T=0:0.5:length(pics);    
end

T=double(T);

 for i=1:length(pics)
IMG(:,:,i)=pics{i};
 end
pics={};

Images=double(IMG);
RawImg=Images(:,:,1);


%% CODE TO MODIFY TCs
st=1;
ed=size(Images,3);
% ed=150;
T=st:0.5:ed;  

% 
images=Images(:,:,st:ed);
DataOut.images=images;
% % % % 

%% DECLARING IMAGE PROPERTIES

sx=size(images,1);
sy=size(images,2);
sz=size(images,3);

%% MASKING IMAGES

if isempty(mask);
    figure
imshow(mat2gray(RawImg));
%numAreas=input('Enter number of areas');
numAreas=1;
mask=zeros(sx,sy);

% level=graythresh(mat2gray(RawImg));
% mask=im2bw(mat2gray(RawImg),level);
disp('select islet')
for i=1:numAreas
bf=imfreehand();
useMask=createMask(bf);
mask=mask+useMask;
end
end
labelbw=bwlabel(mask,4);
numAreas=max(labelbw(:));
currArea=1;
imagesRaw=images;
%%%%%
 h = fspecial('average', 5);
for i=1:sz
%    images(:,:,i)=imfilter(images(:,:,i),h,'replicate').*mask;
   images(:,:,i)=images(:,:,i).*mask;
%    images(:,:,i)=medfilt2(images(:,:,i),[5 5]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Area=size(nonzeros(mask),1);

%% CALCULATING MAPS
timeVec=T;
Fs=mean(diff(timeVec));
Fs=1/Fs;

for i=1:sz
    Mean(i)=mean2(images(:,:,i));
end

oscMap=zeros(sx,sy);
peakMap=zeros(sx,sy);
troughMap=zeros(sx,sy);
timeMat=zeros(sx,sy);
peakAmpMap=zeros(sx,sy);
corrMap=nan(sx,sy);
coeffcorrMap=zeros(sx,sy);
t=1:sz;

while currArea<=numAreas
    
[x,y]=find(mask==currArea);

for i=1:length(x)


               dat=images(x(i)-3:x(i)+3,y(i)-3:y(i)+3,:);

               dat=mean(dat,2);
               dat=mean(dat,1);
               dat=dat(:);
%                size(dat)
%                size(t)
%                f=polyfit(t',dat,2);
%                feval=polyval(f,t');
            
            
%                dat=dat-feval;
           
%                dat=detrend(dat);

               [peakLoc,peakAmp]=peakfinder(dat,(max(dat)-min(dat))/5);
               [trough,troughAmp]=peakfinder(dat,(max(dat)-min(dat))/5,[],-1);
               numpeaks=size(peakLoc,1);
               peakMap(x(i),y(i))=numpeaks;
               troughMap(x(i),y(i))=size(trough,1);
               if length(peakLoc)>1
               timeMat(x(i),y(i))=mean(diff(peakLoc,1));
               peakAmpMap(x(i),y(i))=(mean(peakAmp)-mean(troughAmp));
               end
               if length(peakLoc)==1
                    timeMat(x(i),y(i))=peakLoc;   
                    peakAmpMap(x(i),y(i))=peakAmp-mean(troughAmp);
               end
               if isempty(peakLoc)    
                   timeMat(x(i),y(i))=0;
                   peakAmpMap(x(i),y(i))=0;
               end
               
               
              [c lags]=xcov(dat,Mean,'coeff');
              [maxCV maxCL]=max(c);

               if ~isnan(dat(1))
               corrMap(x(i),y(i))=lags(maxCL)/Fs;
               coeffcorrMap(x(i),y(i))= maxCV;
               end

end
currArea=currArea+1;
end
DataOut.xInfo=length(x);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%ACTIVITY CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% APPLYING "SILENT CELL" THRESHOLD
DataOut.PeakAmpMap=peakAmpMap;

%%%%%%Selecting "silent cell" from peak amp map%%%%%%%%%%
% if selectBackgroundBoolean
    if isempty(silentcell);
    figure;
    imagesc(peakAmpMap)
    disp('select background cell')
    bkgcell=createMask(imfreehand);
    disp('select silent cell')
    silentcell=createMask(imfreehand);
  
    end
    BKGMap=peakAmpMap.*bkgcell;
    bkgthresh=max(max(BKGMap));
    SCMap=peakAmpMap.*silentcell;
    scthresh=max(max(SCMap));

    
    

    DataOut.silentcell=silentcell;
    DataOut.bkgcell=bkgcell;

  
  DataOut.SCthresh=scthresh;
  DataOut.BKGthresh=bkgthresh;

    %%%% using the silent cell threshold to remove areas with NO signal %%%%%%

    [NoSignalMask, NoSignalArea]=RemovingAreaswNoSignal(peakAmpMap,bkgthresh);

% else
%     NoSignalMask = zeros(sx, sy);
%     NoSignalArea = 0;
% 
% end
  


%% CALCULATING ACTIVITY and GENERATING ACTIVITY MAPS

%%%%%% thresholding areas above 2*silent cell threshold%%%%%%%%%%%%%
ActMap=zeros(sx,sy);    
% ActMap(find(mask>0))=1;
ActMap(~mask)=0;
ActMap(NoSignalMask)=0;
ActMap(find(peakAmpMap>=2.0*scthresh))=1; %---- ORIGINAL 
% ActMap(find(peakAmpMap>=1.5*scthresh))=1;
ActMask=logical(ActMap);


% ActMap=bwlabeln(im2bw(ActMap));
% STATS=regionprops(ActMap,'area');




%% MODIFYING ACTIVITY MAPS
 
  colors = [
1 0 1 % First element = purple
0 0 1 % blue
0 1 0 % green
1 1 0 % yellow
1 .65 0 % orange
1 0 0]; % red

n = 256; % size of new color map
m = size(colors,1);
t0 = linspace(0,1,m)';
t = linspace(0,1,n)';
r = interp1(t0,colors(:,1),t);
g = interp1(t0,colors(:,2),t);
b = interp1(t0,colors(:,3),t);
cmap = [r,g,b];

% CoorArea=ActMask*CorrMap;

CoorMap=zeros(sx,sy);
lCoorMap=zeros(sx,sy);

%%%%%% removing regions in the no signal mask, outside of the activity mask, 
%and the image mask%%%%%%%%%%%%
CoorMap=coeffcorrMap;
CorrMap_Final=CoorMap;

CoorMap(~ActMask)=0;




%DataOut.CoorMapOver20_i=CoorMap;

iCoorMap=CoorMap;

lCoorMap=bwlabeln(iCoorMap,4);
%figure;
%imshow(lCoorMap)


%%%%%%%%%%%%removing extremely small areas from maps %%%%%%%%%%%%%%%%%%%%%
STATS=regionprops(lCoorMap,'area');
    stats=0;
for k=1:length(STATS)
      stats(k)=STATS(k).Area; 
   end
   nostats=find(stats<100);
   yesstats=find(stats>100);
%  
% CorrMap=iCoorMap;
for i=1:length(nostats)
CoorMap(find(lCoorMap==nostats(i)))=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CoorMap=CorrMap;
%DataOut.CoorMapOver20_f=CorrMap;

% CoorMap=CorrMap;

%% PLOTTING ACTIVITY MAP


%h = fspecial('average', 2); 
%CoorMap=imfilter(round(CoorMap,1),h,'replicate');
ActiveArea=size(nonzeros(ActMask),1);
% NoNANCorrMap=CoorMap;
CoorMap(find(CoorMap==0))=NaN;
pause(2)
%imoverlayNurin(mat2gray(imagesRaw(:,:,1)),medfilt2(CoorMap,[10, 10], 'symmetric'),[0, 1],[],cmap,0.3)
imoverlayNurin(mat2gray(imagesRaw(:,:,1)),CoorMap,[0, 1],[],cmap,0.3)
colorbar
title('Activity Map')
ActMask=logical(ActMap);
CorrMap(~ActMask)=NaN;


ActMap(~mask)=NaN;
ActMap(NoSignalMask)=NaN;
CorrMap(NoSignalMask)=NaN;

%
% pause(2)
% imoverlay(mat2gray(imagesRaw(:,:,1)),ActMap,[],[],cmap,0.3)
% colorbar

%  h = fspecial('average', 2);
% pause(2)
% imoverlay(mat2gray(imagesRaw(:,:,1)),imfilter(round(CorrMap,1),h,'replicate'),[0, 1],[],cmap,0.3);
% colorbar

% DataOut.round=round(CorrMap,1);

% DataOut.RawImages=imagesRaw;


%% CALCULATING PERCENT ACTIVE AREA

DataOut.RatioActive=ActiveArea/(Area-NoSignalArea);


tdTomMask(~mask)=0;
tdTomMask=imfilter(tdTomMask, [2 2], 'replicate');

tdTomAct=ActMask;
tdTomAct(~tdTomMask)=0;

tdTomCorr=CorrMap_Final;
tdTomCorr(~tdTomMask)=0;
DataOut.CorrCoeffTdPositive=nanmean(nonzeros(tdTomCorr));
DataOut.CorrCoeffActiveAreaTdPositive=nanmean(nonzeros(tdTomCorr(ActMask)));
DataOut.tdPositiveCorrMap=tdTomCorr;
tdTomCorr(~ActMask)=0;
DataOut.tdPositiveCorrMapActiveArea=tdTomCorr;

figure; imagesc(tdTomAct)
DataOut.RatioActiveTdPositive=size(nonzeros(tdTomAct),1)./size(nonzeros(tdTomMask),1);
DataOut.tdPositiveActMap=tdTomAct;
%% SAVING ACTIVITY VARIABLES

DataOut.ActMap=ActMap;
DataOut.ActMapValues = CoorMap;
%DataOut.CorrMap=coeffcorrMap;
DataOut.mask=mask;
DataOut.Location=data.Location;
DataOut.NoSignalMask=NoSignalMask;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *%%%%%%%%%%%%%%%%%%%%%%%%CORRELATION CALCULATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SELECTING THREE CELLS AND SAVING TIME COURSES

sz_2=ed-st+1;

t=1:sz;

if length(fieldnames(data))==1

figure;
imagesc(peakAmpMap)
disp('select cell 1 for correlation analysis')
cell1=createMask(imfreehand);
maskstack = repmat(cell1,[1 1 sz_2]);
maskedStack = images .* maskstack;
cell1TC=shiftdim(mean(mean(maskedStack,2),1),1);
f=polyfit(t,cell1TC,2);
cell1TC_r=cell1TC./polyval(f,t);
cell1TC=medfilt1(cell1TC_r,10);
% 
disp('select cell 2 for correlation analysis')
cell2=createMask(imfreehand);
maskstack = repmat(cell2,[1 1 sz_2]);
maskedStack = images .* maskstack;
cell2TC=shiftdim(mean(mean(maskedStack,2),1),1);
f=polyfit(t,cell2TC,2);
cell2TC_r=cell2TC./polyval(f,t);
cell2TC=medfilt1(cell2TC_r,10);

disp('select cell 3 for correlation analysis')
cell3=createMask(imfreehand);
maskstack = repmat(cell3,[1 1 sz_2]);
maskedStack = images .* maskstack;
cell3TC=shiftdim(mean(mean(maskedStack,2),1),1);
f=polyfit(t,cell3TC,2);
cell3TC_r=cell3TC./polyval(f,t);
cell3TC=medfilt1(cell3TC_r,10);

figure
plot(cell1TC_r)
hold on
plot(cell2TC_r)
plot(cell3TC_r)
legend('cell 1','cell 2','cell 3')
title('timecourses of selected cells')
ylim([0.8 1.3])
hold off

else

   CorrCellTCs=data.CorrCellTCs;
   DataOut.CorrCellTCs_i=CorrCellTCs;
   cell1TC=CorrCellTCs(1,:);
   cell2TC=CorrCellTCs(2,:);
   cell3TC=CorrCellTCs(3,:);
   
end

%% GENERATING CORRELATION MAPS 

% labelbw=bwlabel(mask,4);
% numAreas=max(labelbw(:));
% currArea=1;

clear STATS stats nostats yesstats val loc


[sx,sy,sz]=size(images);

CorrMap1=zeros(sx,sy);
CorrMap2=zeros(sx,sy);
CorrMap3=zeros(sx,sy);

PerOnMap=zeros(sx,sy);

t=1:sz;

intermediateMap = medfilt2(peakAmpMap,[10,10]);
intermediateMap(~mask) = 0;
intermediateMap(NoSignalMask) = 0;

CoorMapBW=(intermediateMap./mean(nonzeros(intermediateMap)));
CoorMapBW(~mask) = NaN;
CoorMapBW(NoSignalMask) = NaN;

LogicalMap=zeros(sx,sy);    
LogicalMap(find(mask==1))=1;
LogicalMap(NoSignalMask)=0;

maskCorr=zeros(sx,sy);
maskCorr(find(LogicalMap==1))=1;
maskCorr=repmat(maskCorr,[1 1 sz]);
Ims=maskCorr.*images;
% 
% dat=shiftdim(mean(mean(Ims(:,:,:),2),1),1);
% f=polyfit(t,dat,2);
% dat=dat./polyval(f,t);
% dat=medfilt1(dat,10);
% 

[x,y]=find(LogicalMap==1);
% [x,y]=find(mask>0);

t=1:sz;
% figure; hold on;
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
    dat=dat./feval;
     
   if ActMask(x(i),y(i))==1
        PerOnMap(x(i),y(i))=DutyCycleCalculator(dat); 
   else
        PerOnMap(x(i), y(i))=0;
   end
    
    dat=medfilt1(dat,10);
    dat=dat';
     [c lags]=xcov(dat,cell1TC, 'coeff');
     [maxCV maxCL]=max(c);
     if ~isnan(dat(1)) 
        % corrMap(x(i),y(i))=lags(maxCL)/Fs;
        CorrMap1(x(i),y(i))= maxCV;
        CorrMap1_2(x(i),y(i))= lags(maxCL);
     end

        [c lags]=xcov(dat,cell2TC, 'coeff');
        [maxCV maxCL]=max(c);
        if ~isnan(dat(1)) 
            CorrMap2(x(i),y(i))= maxCV;
            CorrMap2_2(x(i),y(i))= lags(maxCL);
        end

        [c lags]=xcov(dat,cell3TC, 'coeff');
        [maxCV maxCL]=max(c);
        if ~isnan(dat(1)) 
            CorrMap3(x(i),y(i))= maxCV;
            CorrMap3_2(x(i),y(i))= lags(maxCL);
        end

%         peron=DutyCycleCalculator(dat);
%         PerOn(i)=peron;

end

% dutycycle=mean(PerOn); 
% DataOut.PercentOn=PerOn;
% DataOut.DutyCycle=dutycycle;

CorrMap1=imfilter(round(CorrMap1,1),h,'replicate');
CorrMap2=imfilter(round(CorrMap2,1),h,'replicate');
CorrMap3=imfilter(round(CorrMap3,1),h,'replicate');
CorrMap1(~LogicalMap)=NaN; CorrMap2(~LogicalMap)=NaN; CorrMap3(~LogicalMap)=NaN;

figure('rend','painters','pos',[500 10 400 1000])
subplot(3, 1, 1); imagesc(CorrMap1);colorbar; title('Correlation Map for Cell 1');
subplot(3, 1, 2); imagesc(CorrMap2);colorbar; title('Correlation Map for Cell 2');
subplot(3, 1, 3); imagesc(CorrMap3);colorbar; title('Correlation Map for Cell 3');

DataOut.CorrCellTCs=[cell1TC; cell2TC; cell3TC];
DataOut.CorrCellMaps(:,:,1)=CorrMap1;
DataOut.CorrCellMaps(:,:,2)=CorrMap2;
DataOut.CorrCellMaps(:,:,3)=CorrMap3;

DataOut.CorrCellMaps_2(:,:,1)=CorrMap1_2;
DataOut.CorrCellMaps_2(:,:,2)=CorrMap2_2;
DataOut.CorrCellMaps_2(:,:,3)=CorrMap3_2;


DataOut.PerOnMap=PerOnMap;
DataOut.DutyCycle=mean(nonzeros(PerOnMap));




%% THRESHOLDING DOWN LOW-CORELLATED REGIONS
% 
% CorrMap1(find(CorrMap1<0.1))=0;
 ACorr1=size(nonzeros(CorrMap1.*~isnan(CorrMap1)),1);

% CorrMap2(find(CorrMap2<0.1))=0;
 ACorr2=size(nonzeros(CorrMap2.*~isnan(CorrMap2)),1);
% 
% CorrMap3(find(CorrMap3<0.1))=0;
 ACorr3=size(nonzeros(CorrMap3.*~isnan(CorrMap3)),1);
 
avgCorr1 = nanmean(nanmean(CorrMap1));
%avgCorr1_nozeros = nanmean(nonzeros(CorrMap1));
%sizeCorr1 = size(~isnan(CorrMap1));
avgCorr2 = nanmean(nanmean(CorrMap2));
%avgCorr2_nozeros = nanmean(nonzeros(CorrMap2));
%sizeCorr2 = size(~isnan(CorrMap2));
avgCorr3 = nanmean(nanmean(CorrMap3));
%avgCorr3_nozeros = nanmean(nonzeros(CorrMap3));
%sizeCorr3 = size(~isnan(CorrMap3));

DataOut.AvgCorr = [avgCorr1 avgCorr2 avgCorr3];
%DataOut.AvgCorr_nozeros = [avgCorr1_nozeros avgCorr2_nozeros avgCorr3_nozeros];



%% CALCULATING CORRELATIONS



names=['CorrMap1'; 'CorrMap2'; 'CorrMap3'];
%%%%%%%%%%%%%%%% Plotting Corr Maps %%%%%%%%%%%%%%%%%%%%%%
% pause(2)
% imoverlayNurin(mat2gray(RawImg),CorrMap1,[0,1],[],cmap,0.3)
% colorbar
% title 'Correlation Map for Cell 1'
% pause(2)
% 
% imoverlayNurin(mat2gray(RawImg),CorrMap2,[0,1],[],cmap,0.3)
% colorbar
% title 'Correlation Map for Cell 2'
% pause(2)
% 
% imoverlayNurin(mat2gray(RawImg),CorrMap3,[0,1],[],cmap,0.3)
% colorbar
% title 'Correlation Map for Cell 3'
% pause(2)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
To30_1=zeros(sx,sy);
To50_1=zeros(sx,sy);
To75_1=zeros(sx,sy);

To30_2=zeros(sx,sy);
To50_2=zeros(sx,sy);
To75_2=zeros(sx,sy);

To30_3=zeros(sx,sy);
To50_3=zeros(sx,sy);
To75_3=zeros(sx,sy);

 To30_1(find(CorrMap1<0.3 & CorrMap1>0))=0.25;
 sTo30(1)=size(nonzeros(To30_1),1);
 To50_1(find(CorrMap1<0.7 & CorrMap1>=0.3))=0.5;
 sTo50(1)=size(nonzeros(To50_1),1);
 To70_1(find(CorrMap1<=1.0 & CorrMap1>=0.7))=1;
 sTo70(1)=size(nonzeros(To70_1),1);

 To30_2(find(CorrMap2<0.3 & CorrMap2>0))=0.25;
 sTo30(2)=size(nonzeros(To30_2),1);
 To50_2(find(CorrMap2<0.7 & CorrMap2>=0.3))=0.5;
 sTo50(2)=size(nonzeros(To50_2),1);
 To70_2(find(CorrMap2<=1.0 & CorrMap2>=0.7))=1;
 sTo70(2)=size(nonzeros(To70_2),1);
 
 To30_3(find(CorrMap3<0.3 & CorrMap3>0))=0.25;
 sTo30(3)=size(nonzeros(To30_3),1);
 To50_3(find(CorrMap3<0.7 & CorrMap3>=0.3))=0.5;
 sTo50(3)=size(nonzeros(To50_3),1);
 To70_3(find(CorrMap3<=1.0 & CorrMap3>=0.7))=1;
 sTo70(3)=size(nonzeros(To70_3),1);

 

ACorrs=[sTo30; sTo50; sTo70];
ACorrs=ACorrs/(Area-NoSignalArea);
DataOut.ACorrs=ACorrs;

%[ACorr,ind]=sort([sTo30; sTo50; sTo70],'descend');
RatioCorr_1=(sTo70(1))./(Area-NoSignalArea);
RatioCorr_2=(sTo70(2))./(Area-NoSignalArea);
RatioCorr_3=(sTo70(3))./(Area-NoSignalArea);

CorrCoeffs=[RatioCorr_1 RatioCorr_2 RatioCorr_3];
DataOut.CorrCoeffs=CorrCoeffs;

sCC=sort(CorrCoeffs, 'descend');
HiCorrCoeff=sCC(1);
DataOut.HiCorrCoeff=HiCorrCoeff;

disp('End of analysis.')

end