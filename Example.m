

oLocation=pwd;
Location='C:\Users\Undergrads\Desktop\Nurin\Videos_test2\';
cd(Location)

%%
files = dir('*.czi');
for i=1:length(files)
    
  [pathname,filename,extension] = fileparts(files(i).name);
  % the new name, e.g.
  [sx sy]=size(filename);
  newFilename = filename(2:end);
  % rename the file
  movefile([filename extension], [newFilename extension]);
end

%%

filePattern = fullfile(Location)
list = dir (filePattern );
[B ,index] = sortrows({list.name});
 
[sx,sy]=size(B);
B=B(3:sy-1);

filePattern_f = fullfile(filePattern, B(:));


for i=1:length(B)
ff=char(filePattern_f(i));
[pathstr,ne{i},ext] = fileparts(ff);
end


cd(oLocation)
for i=1:length(ne)
Loc=char(filePattern_f(i));
nam=char(ne(i));
Analysis.(char(ne{i}))=Run(Loc);
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
% WT_11G_3MH(4,:)=[WT_11G_3MH_I1.RatioActive WT_11G_3MH_I1.HiCorrCoeff WT_11G_3MH_I1.DutyCycle];
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

    