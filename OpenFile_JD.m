function [data] = OpenFile_JD( data )
% open file with bfopen and remove any unneeded info
if data.newmask == 1
    data.mask=[];
end
R=bfopen(data.Location);
pics=R{1};
pics=pics(:,1);% remove directory info

%get time slices
try R{4}.getPlaneDeltaT(0,2); 
    if ~isempty(R{4}.getPlaneDeltaT(0,0))
        for i=1:length(pics)
            T(i)=R{4}.getPlaneDeltaT(0,i-1);
        end
    else
        T=1:size(pics,3);
    end
catch
    % for i=1:length(pics)
    % TSize(i)=R{4}.getPlaneDeltaT(i-1,0);
    % end    
    T=0:0.5:length(pics);    
end

T=double(T);

%change variable for pics
for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
pics={};
data.images = IMG;
data.TSize = T;
end

