% dir2='d:\Science\UUtrecht\Experimental\22_10_18_contact_sites_Life_take2\Processed\ROI_3c\overlay_files\';%
dir2='d:\Science\UUtrecht\paper_with_gareth\Figures_systematic_attempt\original_data_sample\3D_FIB_CLSM_volumes\case_1_EDGE_sample_standard\Mlb_process\';
file_template='C1-Stack_of_overlays_2channels';%C4-overlay_tiff C3-overlay_1_tiff
start=1;N=7;%starting frame, number of frames
oldFolder=cd(dir2);

% % number=SnV_n2s(0+start);
number=SnV_n2fixed_length_string(0+start,2);

image=imread([file_template,number,'.tif']);
[h0,w0] = size(image);
Stack=uint8(zeros(h0,w0,N));
Stack(:,:,1)=image;
h=waitbar(0,'Loading images...');
for j=(1:(N-1)) 
% %     number=SnV_n2s(j+start);
    number=SnV_n2fixed_length_string(j+start,2);
% %     image=imread([file_template1,number,'.tif']);
    Stack(:,:,j+1)=imread([file_template,number,'.tif']);    
    waitbar(j/(N-1))
end
close(h);
cd(oldFolder);