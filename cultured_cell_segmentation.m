addpath 'Z:\NoRI\NoRI_code\nori.machine.code\'
addpath 'N:\NoRI\NoRI_code\srs.code\SRS functions\'
addpath(genpath('M:\MATLAB\xiliFunctions'))
cd 'N:\NoRI\Seungeun\20231110 cell segmentation code\'

%% 
clear all
basepath=pwd;
imname='nori_image_HeLa.tif';
im= fasttifread(strcat(basepath, filesep, imname));

coef = 8192; % NoRI data scaling factor 
lipid_density = 1.010101; % g/ml,
protein_density = 1.364256; % g/ml,
dx = 212.132/512; % micrometer/pixel

% Parse image channels
protein = double(im(:,:,1))*protein_density/coef;
lipid = double(im(:,:,2))*lipid_density/coef;
figure(1),imagesc(protein,[0 0.1]),axis image
figure(2),imagesc(lipid,[0 0.1]),axis image

% Threshold values 
cellthres=400*protein_density/coef;          % Threshold for protein density of cells
nuthres=100*protein_density/coef;            % Threshold for protein density of nucleus
cell_size=500;          % Size threshold for cell
nucleus_size=250;       % Size threshold for nucleus
nucleoli_size=3;        % Size threshold for nucleoli

% Cell body mask
bw = protein>cellthres; % 
bw = imclose(bw,strel('disk',4));
bw = imopen(bw,strel('disk',4));
bw = bwareaopen(bw,500);
bw_cell = bw;
figure(1),imagesc(labeloverlay(rescale(protein,0, 4),imdilate(bwperim(bw),strel('disk',3)),'colormap','autumn')), axis image

% Nucleus mask
bw = lipid<nuthres; 
bw = imclose(bw,strel('disk',2));
bw = bw.*bw_cell;
bw = imopen(bw,strel('disk',4));
bw = imfill(bw,'holes');
bw = bwareaopen(bw,nucleus_size);
bw_large = bwareaopen(bw,5000);
bw = bw - bw_large;
bw_nu = bw;
figure(2),imagesc(labeloverlay(rescale(medfilt2(lipid,[3 3]),0, 5),imdilate(bwperim(bw_nu),strel('disk',2)),'colormap','autumn')),axis image


% Segment cell body with nucleus mask
I1 = medfilt2(protein,[2 2]);
I1 = imopen(I1,strel('disk',10));
D = imimposemin(-I1, bw_nu);
L = watershed(D);
L = L&bw_cell;
L = imfill(L,'holes');
L = bwareaopen(L,cell_size);
L = bwlabel(L);
L_nu = bw_nu.*L;
L_cell = L;
figure(3),imagesc(labeloverlay(rescale(medfilt2(lipid,[3 3]),0, 10),L_cell,'colormap','autumn')), axis image
figure(4),imagesc(labeloverlay(rescale(medfilt2(lipid,[3 3]),0, 10),L_nu,'colormap','autumn')), axis image

% Nucleoli mask (simple version for homogeneous population)
protein_med = medfilt2(protein,[3 3]);
thresh_ncl = multithresh( protein_med(find(bw_nu)) ,2);
bw_ncl=(protein>thresh_ncl(end)).*bw_nu; 
bw_ncl = bwareaopen(bw_ncl,nucleoli_size);
figure(2),imagesc(labeloverlay(rescale(medfilt2(protein,[3 3]),0, 2),bwperim(bw_ncl),'colormap','autumn')),axis image

% Nucleoli mask (auto-threshold applied to each nucleus)
bw_ncl=zeros(size(protein));
for k=1:max(L_nu(:))
    if nnz(L_nu==k)>0
        nu_property = regionprops(uint16(L_nu==k),protein,'PixelIdxList','Solidity');
        thresh_ncl = multithresh( protein([nu_property.PixelIdxList]) ,2);
        if nu_property.Solidity>0.7
            bw_ncl_temp=(protein>thresh_ncl(end)).*(L_nu==k);
            bw_ncl_temp = bwareaopen(bw_ncl_temp,nucleoli_size);
            bw_ncl=bw_ncl+bw_ncl_temp;
        end
    end
end
figure(2),imagesc(labeloverlay(rescale(protein_med,0, 2),bwperim(bw_ncl),'colormap','autumn')), axis image

% Nucleoli and Nucleoplasm segmentation
L_ncl=L_nu.*bw_ncl;
L_nplasm=L_nu.*(1-bw_ncl);

% Export segmentation
segpath = strcat(basepath,filesep,'Output\');
mkdir(segpath)
save(strcat(segpath, imname(1:end-4),'_seg.mat'),'L_cell','L_nu','L_ncl','L_nplasm')%
writeTIFF(cat(3,protein,lipid,single(L_cell>0),single(L_nu>0),single(L_ncl>0)), strcat(segpath, imname(1:end-4),'_seg.tif'),'w');