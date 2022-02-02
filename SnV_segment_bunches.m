% Manual SnV segmentation using freehand selection, at one-in-ten slices
% dim=3;% dimension of slicing
bunching=20;%N=500;
% i=1;%counter

order=[1,3,2];
Stack=permute(Stack_0,order);%.*uint8(~(Stack_V_OR))|Stack_ER_321_1|Stack_Mito_OR|Stack_ER_312_1))
Stack_segmented=zeros(size(Stack),'logical');
N=size(Stack,3);
start=1;
close all;figure;
for i=(1:(ceil((N-start)/bunching)))%size(Stack,dim)
    I_bunch=Stack(:,:,(start+(i-1)*bunching):min(start+i*bunching,N))...
    .*uint8(~Stack_segmented(:,:,(start+(i-1)*bunching):min(start+i*bunching,N)));
    %% Image for user
    I_0=mean(I_bunch,3);%becomes double!
    cla;imagesc(I_0);%%-double(min(I_bunch,[],3))
    axis equal;axis tight;title([num2str(i*bunching/N*100),' %']);
    colormap gray;%caxis([0 100]);
    %% Segmentation by user
    mask_bunch=zeros(size(I_bunch),'logical');I_detected=mask_bunch;I_t=mask_bunch;%I_manual=mask;
    size_check=1; % variable for ending process
    while size_check    %% select all the membrane pieces, untill they are over
    ans=[];
     while isempty(ans) %% should allow to redraw
      h = imfreehand();wait(h)% drawfreehand(gca,'Multiclick',true);        
     end
    I_manual=createMask(h); %% user selected area
     if (sum(double(I_manual(:)))>10)% click twice to go for the next bunch
      treshold=multithresh(I_0.*double(I_manual),2); %vesicle 1, ER - 2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      %% Auto-process           
       parfor j=1:bunching
        I_t=I_bunch(:,:,j).*uint8((I_bunch(:,:,j)>treshold(2)).*I_manual);%tresholded vesicle 1, ER - 2 
        I_detected=imerode(imfill( edge( I_t,'Canny',0.7,5)+(I_t>0),'holes'),ones(3));% 
        mask_bunch(:,:,j)=mask_bunch(:,:,j)+I_detected;
       end
      imagesc(double(~max(mask_bunch,[],3)).*I_0);%updated image %-double(min(I_bunch,[],3))
      axis equal;axis tight;%ready for new selection
%     caxis([0 100])
      title([num2str((i+1)*bunching/N*100),' %']);
     else
      size_check=0;
     end    
    end
    Stack_segmented(:,:,(start+(i-1)*bunching):(start+i*bunching))=mask_bunch;
end
Stack_segmented=ipermute(Stack_segmented,order);
%% Visualisations:
% % imagesc(max(Stack_ER_OR,[],3));axis equal;axis tight;
figure;

hold on;grid on;
p = patch(isosurface(double(smooth3(fliplr(Stack_ER_123(:,:,378:-1:1))))));
p.FaceColor = 'b';
p.EdgeColor = 'none';
daspect([1 1 27/128]);
axis equal;

% lights=camlight;% important for recording the videos
% lighting phong

%%%% axis([1 1024 1 827 1 20])
%% Save results
h = waitbar(0,'Saving...');
% % %%%%%%%Stack2=255*(Stack2-min(Stack2(:)))/(max(Stack2(:))-min(Stack2(:)));
% Stack2=uint8(smooth3(Stack_segmented));
file2='MLb_ER_XY.tiff';
delete(file2)
for i=1:size(Stack_ER_123,3)
    imwrite(uint8(Stack_ER_123(:,:,i)),file2,'WriteMode','append')
    waitbar(i/size(Stack_ER_123,3))
end
close(h)