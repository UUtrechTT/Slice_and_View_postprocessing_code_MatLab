%% filtering and finding particles in 2D/3D version
N2=700;
start2=0;
%% Particle templates
dim_p=[38,38,10];
Prtcl0=zeros(dim_p);

[X_p,Y_p,Z_p]=meshgrid(2*((0:dim_p(1)-1)/(dim_p(1)-1)-0.5),2*((0:dim_p(2)-1)/(dim_p(2)-1)-0.5),2*((0:dim_p(3)-1)/(dim_p(3)-1)-0.5));%coordinate grid from -1 to 1
core=0.03;%0.09 at 56
Core0=Prtcl0+...
    1*double(((X_p).^2+(Y_p).^2+(Z_p).^2)<core).*( ( -((X_p).^2+(Y_p).^2+(Z_p).^2)+core)/core).^0.5;
Prtcl0=Prtcl0+...
    1*double(((X_p).^2+(Y_p).^2+(Z_p).^2)<1.0)... %shell
    +0*Core0... %core
    +0*double(((X_p).^2+(Y_p).^2+(Z_p).^2)>0.7).*double(((X_p).^2+(Y_p).^2+(Z_p).^2)<1.1);% rim

Prtcl=squeeze(Prtcl0(:,:,ceil(dim_p(3)/2)));%2D version
Core=squeeze(Core0(:,:,ceil(dim_p(3)/2)));%2D version
% Prtcl=Prtcl-ordfilt2(Prtcl,1,true(2),'symmetric');

Core=edge(Core,'canny',0.1,3);
Core=padarray(Core,[(h0-dim_p(1))/2 (w0-dim_p(2))/2], 0);
Core_S=(2*pi)^4.*conj(fft2(Core))./sum(abs(Core(:)).^2).^0.5;

% Prtcl=edge(Prtcl,'canny',0.1,3);
Prtcl=padarray(Prtcl,[(h0-dim_p(1))/2 (w0-dim_p(2))/2], 0);%size of frame
Prtcl_S=fft2(Prtcl);

dim_new=size(Prtcl);%in case of small size template
%% Reisz transform in k space
%
freq_x=repmat((1:dim_new(1))/(dim_new(1)+2)-0.5,dim_new(2),1)';%frequency space
freq_y=repmat((1:dim_new(2))/(dim_new(2)+2)'-0.5,dim_new(1),1);%

% % kernel_x=(fftshift(-1i*freq_x.*(freq_x.^2+freq_y.^2).^-0.5));kernel_x(isnan(kernel_x))=0;
% % kernel_y=(fftshift(-1i*freq_y.*(freq_x.^2+freq_y.^2).^-0.5));kernel_y(isnan(kernel_y))=0;
% % kernel_xx=kernel_x.^2;kernel_yy=kernel_y.^2;%imagesc((abs(fftshift(kernel_xx))));axis([570 584 1451 1465])
% Prtcl_S=Prtcl_S.*fftshift( 1-0.01*(0.01+freq_x.^2+freq_y.^2).^-1 );%figure;imagesc(fftshift(abs(Prtcl_S)))
% % Prtcl_xx=-(Prtcl_S.*kernel_xx);
% % Prtcl_yy=-(Prtcl_S.*kernel_yy);
% % Prtcl_xy=-(Prtcl_S.*kernel_x.*kernel_y);
% imagesc(abs(Prtcl_S))
%%
% figure;
Stack2=zeros(h0,w0,N2);
h = waitbar(0,'Searching...');
tic;
treshold_core=220;
core_area_treshold=[100 1000];%[10 75]
treshold_shell_1=50;
treshold_shell_2=150;
shell_area=[50 150];

treshold_bckg=80;
particle2D_size = strel('disk',dim_p(1)/2,4);
% se = strel('disk',11);
n=zeros(1,N2);Ac=n;
% % figure;
for i=1:N2
%%
    image=Stack(:,:,i+start2);% .*uint8(~bwareafilt((Stack_T2(:,:,i+start2)>130),1))
    
    cores=bwareafilt(image>treshold_core,core_area_treshold);
% %cores=(bwperim(cores,8));%%perimeter->inversion->size selection; filters out candidates coinsiding with their perimeters
% %cores=bwareafilt(not(cores),[1 core_area_treshold(2)]);%%.*imregionalmax(image,8)
    
    s=regionprops(cores, 'Area','centroid');% just to make sure we have the same numeration
    core_areas=  ceil(cat(1, s.Area    ));
    core_centers=ceil(cat(1, s.Centroid));
% %cores=cores.*imregionalmax(image,8);
% %cores=bwpropfilt(cores,'ConvexArea',core_area,4);
    ind=cor2ind2D(core_centers,[w0 h0]);
    ind2=find(image(ind)>treshold_core);%leave only those which have bright centers
    
    core_centers=core_centers(ind2,:);
    core_areas=core_areas(ind2);
    n(i)=numel(ind2);%count of the first picks
    if n(i)>0
        ind=find( (core_centers(:,1)>dim_p(1)/2)&(core_centers(:,1)<h0-dim_p(1)/2)...
                 &(core_centers(:,2)>dim_p(2)/2)&(core_centers(:,2)<w0-dim_p(2)/2) );% not too close to the edge
        core_centers=core_centers(ind,:);
        core_areas=core_areas(ind);
    end
    
    around_cores=double(image).*imdilate(cores,particle2D_size);
% %cores=uint8(cores).*image;    
    number_of_cores=size(core_centers,1);
    n(i)=number_of_cores;
    
% % % %     while number_of_cores>0
% % % %        [~,ind_first_core]=max(core_areas);
% % % %        x1=core_centers(ind_first_core,1);
% % % %        y1=core_centers(ind_first_core,2);
% % % %        
% % % %        local_box=image((x1-dim_p(1)/2):(x1+dim_p(1)/2),(y1-dim_p(2)/2):(x1+dim_p(2)/2));
% % % %        
% % % %        ind2=find((abs(core_centers(:,1)-x1)>dim_p(1)/2)...
% % % %                 &(abs(core_centers(:,2)-y1)>dim_p(2)/2));
% % % %        core_centers=core_centers(ind2,:);%eliminate nearby candidates    
% % % %     end
        
    %cores=bwulterode(image>treshold_core,'cityblock',4)&cores;% detect particle cores 
   
    %shells=bwareafilt(and((around_cores>treshold_shell_1),(around_cores<treshold_shell_2)),shell_area);
    
%     Ac(i)=numel(find(around_cores>0));
% %     S2=fft2(double(shells));% bwpropfilt(,'Eccentricity',[0.0 0.7]) 
% %         corr=fftshift(...
% %      real(ifft2(S2.*Prtcl_S))... %.*sum(abs(S2(:)).^2).^-0.5  .*real(ifft2(Prtcl_S.*S1)).*
% %         );%.*image

    Stack2(:,:,i)=around_cores;%imdilate(cores,se);%.*image
    
    waitbar(i/N2)
end
% Stack2=imdilate(Stack2,strel('ball',dim_p(1)/2,dim_p(3)/2));
close(h)
%%
figure;
subplot(1,2,1);imagesc(squeeze(mean(Stack2,2)));axis equal;axis tight;
subplot(1,2,2);imagesc(+squeeze(max(Stack2(:,:,:),[],3)));axis equal;axis tight;
%%

h = waitbar(0,'Saving...');
Stack2=255*(Stack2-min(Stack2(:)))/(max(Stack2(:))-min(Stack2(:)));
Stack2=uint8(Stack2);
file2='MLb_corr_T1.tiff';
delete(file2)
for i=1:N2
    imwrite(Stack2(:,:,i),file2,'WriteMode','append')
    waitbar(i/N2)
end
close(h)

toc;