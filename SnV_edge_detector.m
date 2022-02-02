%% SnV edge detection

dir0='d:\Science\UUtrecht\Experimental\22_10_18_contact_sites_Life_take2\Processed\ROI_1\';

dir =[dir0,'overlay_files\'];%d:\Experimental\CLEM\05_02_18_HeLa_small_ROI_3D\cell5\cell5_SnV2\TiffImages\';
cd(dir);
dir2=[dir0,'overlay_edges\'];%save here

file_template1='C3-overlay_1_tiff';%

start=0;N=527;%starting frame, number of frames

save=1;
if save
    [~,~,~]=mkdir(dir2);
end
% acquire=false;
filter_images=false;

h = waitbar(0,'Detecting...');

for i=1:size(Stack,2) 
    number=SnV_n2fixed_length_string((i-1)+start,4);
%     image=imread([file_template1,number,'.tif']);%use the first channel for the ROI detection
    image=squeeze(Stack(:,i,:));%
    image_edges=uint8(edge(~bwareafilt(image>50,1).*(image>60),'canny',0.2,5));
    if save
        imwrite(image_edges.*image,[dir2,'edge_',number,'.tif']);
    end
    waitbar(i/(size(Stack,2) -1))
end
close(h);