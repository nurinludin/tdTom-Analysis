function [dat]=Run(Location)
%%
R=bfopen(Location);

pic=R{1};
pic=pic(:,1);

for i=1:length(pic)
    img(:,:,i)=pic{i};
end
pic={};
img1=img;
img=double(img);
% mask=data.mask;

dat.Location=Location;
dat.Img=img(:,:,:);

[MSize, NSize,SSize]=size(img);
count=1;
figure;
for i=1:2:(SSize-1)
    clear BW img1
img1=mat2gray(img(:,:,i));
BW = imbinarize(img1, 'global');
imshowpair(img1,BW,'montage')
areaTom(count)=size(nonzeros(BW),1);

BWStacks(:,:,count)=BW;

count=count+1;

end

count=1;
figure
for i=2:2:SSize

imshow(mat2gray((img(:,:,i))));
matrixMaskLive(:,:,count)=createMask(imfreehand);
% imshow(matrixMaskLive(:,:,i))
areaLive(count)=size(nonzeros(matrixMaskLive(:,:,count)),1);
count=count+1;
end


    ratioMx=areaTom./areaLive;


dat.RatioMxs=mean(ratioMx);
dat.ratioMxs=ratioMx;
dat.areaLive=areaLive;
dat.areaTom=areaTom;
dat.LiveMasks=matrixMaskLive;
dat.tdTomMasks=BWStacks;




end