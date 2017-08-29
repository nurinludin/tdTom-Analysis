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
for i=1:length(ne)
Loc=char(filePattern_f(i));
nam=char(ne(i));
disp(sprintf('Running analysis for image number %d',i)); 
Analysis.(char(ne{i}))=Run(Loc); %,Analysis.(char(ne{i})) );
end




%% Organizing To Matrices
Analysis.Names=(ne);
I=['I1'; 'I2' ;'I3' ;'I4' ;'I5'; 'I6'];

count=1;
for i=1:length(ne)

            data(i,:)=[Analysis.Names(i) Analysis.(char(ne{i})).RatioMxs];

end

Analysis.DataOverview=data;

%% Saving

save('I:\Nurin\JoshData\Analysis\MatlabFiles\042617.mat','Analysis_042617', '-v7.3')

%% Removing a field from the structure
for i=1:length(ne)
    
    if isfield(Analysis.(char(ne{i}))(:),'I3')==1
        for j=1:3
            Analysis.(char(ne{i})).(I(j,:))=rmfield(Analysis.(char(ne{i})).(I(j,:)),'images');
        end
        
    elseif isfield(Analysis.(char(ne{i}))(:),'I2')==1
                  for j=1:2
           Analysis.(char(ne{i})).(I(j,:))=rmfield(Analysis.(char(ne{i})).(I(j,:)),'images');
                  end
%          
    elseif isfield(Analysis.(char(ne{i}))(:),'I1')==1
        j=1;
            Analysis.(char(ne{i})).(I(j,:))=rmfield(Analysis.(char(ne{i})).(I(j,:)),'images');
%        
    end

    
end




%%
% 
% close all
% 
% WT_2G_I1.Location='L:\Images\MH\04-19-2017\WT_2G_I1.czi';
% WT_2G_I1=Run(WT_2G_I1);
% 
% WT_11G_I1.Location='L:\Images\MH\04-19-2017\WT_11G_I1.czi';
% WT_11G_I1=Run(WT_11G_I1);
% 
% WT_11G_3MH_I1.Location='L:\Images\MH\04-19-2017\WT_11G_3MH_I1.czi';
% WT_11G_3MH_I1=Run(WT_11G_3MH_I1);
% 
% WT_2G(4,:)=[WT_2G_I1.RatioActive WT_2G_I1.HiCorrCoeff WT_2G_I1.DutyCycle];
% WT_11G(4,:)=[WT_11G_I1.RatioActive WT_11G_I1.HiCorrCoeff WT_11G_I1.DutyCycle];
% WT_11G_3MH(4,:)=[WT_11G_3MH_I1.RatioActi ve WT_11G_3MH_I1.HiCorrCoeff WT_11G_3MH_I1.DutyCycle];
%   
% WT_data=[WT_2G(:,1) WT_11G(:,1) WT_11G_3MH(:,1) ...
%                   WT_2G(:,2) WT_11G(:,2) WT_11G_3MH(:,2) ...
% WT_2G(:,3) WT_11G(:,3) WT_11G_3MH(:,3)]

%%


%%
% cd(home)
% save('Exp6b_B6.mat', 'B6_*')
% save('Exp6b_WT.mat', 'WT_*')
% save('Exp6b_WT.mat', 'WT_*')

    