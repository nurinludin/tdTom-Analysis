function Analysis=MasterRun(Location)

oLocation=pwd;

cd(Location)

%% Rename Files by character size
% files = dir('*.czi');
% for i=1:length(files)
%     
%   [pathname,filename,extension] = fileparts(files(i).name);
%   % the new name, e.g.
%   [sx sy]=size(filename);
%   newFilename = filename(2:end);
%   % rename the file
%   movefile([filename extension], [newFilename extension]);
% end

%% Retrieving file names



filePattern = fullfile(Location);
list = dir (filePattern );
[B ,index] = sortrows({list.name});
 
[sx,sy]=size(B);
B=B(4:sy-1);

filePattern_f = fullfile(filePattern, B(:));


for i=1:length(B)
ff=char(filePattern_f(i));
[pathstr,ne{i},ext] = fileparts(ff);
end


cd(oLocation)
close all;

ne = natsort(ne);
filePattern_f=natsort(filePattern_f);

%% Running analysis through list
for i=2:6   %:length(ne)
% Loc=char(filePattern_f(i));
disp(sprintf('Running analysis for image number %d',i)); 
nam=char(ne(i))
Analysis.(char(ne{i})).Location=char(filePattern_f(i));
% Analysis.(char(ne{i}))=Run(Analysis_Exp8a.(char(ne{i})));
Analysis_Exp8a.(char(ne{i}))=RunwTdTom(Analysis_Exp8a.(char(ne{i})),Analysis_tdTom.(char(ne{i})).tdTomMasks); %,Analysis.(char(ne{i})) );
end

%%

for i=1:6
    
Analysis_Exp8a.tdTom.Names(i,1)=(ne(i));
Analysis_Exp8a.tdTom.Overview(i,1)=Analysis_tdTom.(char(ne{i})).RatioMxs;
Analysis_Exp8a.tdTom.Overview(i,2)=Analysis_Exp8a.(char(ne{i})).RatioActiveTdPositive;
Analysis_Exp8a.tdTom.Overview(i,3)=Analysis_Exp8a.(char(ne{i})).CorrCoeffTdPositive;

end
%% Organizing Data to arrays

for i=1:length(ne)
    
Analysis.Names(i,1)=(ne(i));
Analysis.Overview(i,1)=Analysis.(char(ne{i})).RatioActive;
Analysis.Overview(i,2)=Analysis.(char(ne{i})).HiCorrCoeff;
Analysis.Overview(i,3)=max(Analysis.(char(ne{i})).AvgCorr);
Analysis.Overview(i,4)=Analysis.(char(ne{i})).DutyCycle;

end


