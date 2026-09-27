%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------
clear all
close all
clc

addpath tools


file_path = './';

load ./data_pack

ESP_eff=3.6000e-04;
N=[200,200];
Ry=4;
tol=0.01
num_dwi=1;
num_chan=52;
num_shot=2;

delta_ky_ini=[0,2; 0,2; 0,2; 3,0; 0,2];
%delta_ky_ini=[0,2; 1,1; 2,0; 3,3; 0,2];

slc_sel=[15];
num_slice=length(slc_sel);

%%
for ii_te=1:5
for ii_rf=1:5
    if 1%subs_mat(ii_rf,ii_te)
    
shot_start = 1+(ii_rf-1)*2;
Del_ky = delta_ky_ini(ii_te,:);
A = fftc(eye(N(1)),1);          
kline_ap=(1+Del_ky(1):Ry:200);
kline_pa=(1+Del_ky(2):Ry:200);
A_ap = A(kline_ap,:);
A_pa = A(kline_pa,:);

PE_line=50;

K_epi_ap=K_epi_ap_all(:,:,:,:,ii_rf,ii_te);
K_epi_pa=K_epi_pa_all(:,:,:,:,ii_rf,ii_te);

sgnl_ap=ifftc(K_epi_ap,1);
sgnl_pa=ifftc(K_epi_pa,1);




  
%
%--------------------------------------------------------------------------
%% MUSSELS: gradient descent with B0 model 
%--------------------------------------------------------------------------
if 1
t_axis_ap = [0:Ry:PE_line*Ry-1] * ESP_eff/2;

t_axis_pa = t_axis_ap(end:-1:1);


winSize = [1,1] * 10;
N(2)=200;
%lambda_msl = 1.25;
%keep = 1:floor(lambda_msl*prod(winSize));
soft_th=0.05;
lambda_msl=0.6;


step_size = 1.5;
num_iter = 200;




img_msl = zeross([N,num_shot,num_slice,num_dwi]);


tic
for nslc = 1% 1:num_slice
    disp(['slice: ', num2str(nslc)])
   % diff = exp(-1i* angle(img_sense_topup(:,:,1,nslc) ./ img_sense_topup(:,:,2,nslc)));
  
    
    % create and store encoding matrices
    AWC_ap = zeross([num_chan*PE_line, N(2), N(1)]);
    AWC_pa = zeross([num_chan*PE_line, N(2), N(1)]);

    AWC_apH = zeross([N(2), PE_line*num_chan, N(1)]);
    AWC_paH = zeross([N(2), PE_line*num_chan, N(1)]);

    for xn = 1:N(1)
        b0 = img_fieldmap(xn,:,nslc) * 2 * pi;

        W_ap = exp(1i * t_axis_ap.' * b0);
        AW_ap = A_ap .* W_ap;      

        W_pa = exp(1i * t_axis_pa.' * b0);
        AW_pa = A_pa .* W_pa;      

        for c = 1:num_chan       
            AWC_ap(1 + (c-1)*PE_line : c*PE_line, :, xn) = AW_ap * diag( sens(xn,:,nslc,c) );
            AWC_pa(1 + (c-1)*PE_line : c*PE_line, :, xn) = AW_pa * diag( sens(xn,:,nslc,c) ); %.* diff(xn,:) );
        end

        AWC_apH(:,:,xn) = AWC_ap(:,:,xn)';
        AWC_paH(:,:,xn) = AWC_pa(:,:,xn)';
    end    
    
   
    for ndwi = 1:num_dwi
        % initialize with sense solution: [helps slightly]
        %im_rec = permute(cat(3, img_sense_ap(:,:,nslc,ndwi), img_sense_pa(:,:,nslc,ndwi)), [2,3,1]);
        im_rec = permute( zeros([N,num_shot]),[2,3,1]);
%       im_rec = permute( squeeze(img_sense(:,:,:,nslc)),[2,3,1]);

        for iter = 1:num_iter
            im_prev = im_rec;
            
            for xn = 1:N(1)
                for ii=1:num_shot
               

                 if ii<=num_shot/2
                rhs_ap = sgnl_ap(xn, :, :, nslc, ii);   
                rhs_ap = rhs_ap(:);  
                im_rec(:,ii,xn) = im_rec(:,ii,xn) - step_size * AWC_apH(:,:,xn) * ( AWC_ap(:,:,xn) * im_rec(:,ii,xn) - rhs_ap );
                 else
                rhs_pa = sgnl_pa(xn, :, :, nslc, ii-num_shot/2);   
                rhs_pa = rhs_pa(:); 
                 im_rec(:,ii,xn) = im_rec(:,ii,xn) - step_size * AWC_paH(:,:,xn) * ( AWC_pa(:,:,xn) * im_rec(:,ii,xn) - rhs_pa );
                 end
                end
            end

 %% low rank           
            im_rec = permute(im_rec, [3,1,2]);

            A = Im2row( fft2call(im_rec), winSize );
            [U, S, V] = svd(A, 'econ');
            S_save(:,iter)=diag(S);
% soft threshold
        %   temp_a=diag(S);
        %    temp_b=find(temp_a(:)<soft_th*temp_a(1));
        %   keep=1:temp_b(1);
        %   keep_count(iter)=keep(end);
% hard threshold
      keep=1:prod(winSize)*num_shot*lambda_msl;
%
            A = U(:,keep) * S(keep,keep) * V(:,keep)';

            k_pocs = Row2im(A, [N, num_shot], winSize);
            im_rec_disp = ifft2call(k_pocs);


%%
            im_rec = permute(im_rec_disp, [2,3,1]);   

     update = rmse(im_prev,im_rec)
    if update < tol
        break
    end
       %     mosaic(im_rec_disp, 1, 2, 6, ['iter: ', num2str(iter), '  update: ', num2str(rmse(im_prev,im_rec))], genCaxis(im_rec), 90)
        end

        img_msl(:,:,:,nslc,ndwi) = permute(im_rec, [3,1,2]);
    end
    
end
toc


%mosaic(mean(abs(img_msl(:,:,:,nslc,1)),3), 1, 1, 6, 'buda', [0,.7e-3], 90)
end

    
%% final
img_final(:,:,:,ii_rf,ii_te)=squeeze(mean(abs(img_msl),3));
close all
clc
ii_rf
ii_te
    end
end
end
%
img_mix=img_final;
%%
%clearvars -except img_final img_final_1 img_final_2 img_final_3 subs_mat img_msl
%%

clear img_recon

 close all
% 1.50 x 1.50 x 0.86
%%
addpath ('library')
addpath ('imagine')
disp('Loading data ....');

Voltfac=1.0;
B1_index=[0.5:0.05:1.5];
B1_index=B1_index*Voltfac;
nrf_d=5;
v_num=5;

B1map_de=1*ones(200,200);

if 0
disp('Load slice profile')
    [s_name, s_path]=uigetfile('*.mat','select slice profile');
    s_raw=strcat(s_path,s_name); 
    load (s_raw)
end

%load ./RF_5x_90thick_1p00t_4p3mm_Verse1_CLv1_B1_Mzsim_scale1p0_MOut    
load ./RF_5x_90thick_1p00t_5p0mm_Verse1_CLv1_B1_Mzsim_T1_1000_TR2p1_MOut
%% load img & V    

img=img_mix;
dic_path='./i_shuffling_dic11_withDic';
load (dic_path)

n_te=size(img,5);

    if 1
  img_ini=reshape(img,[size(img,1),size(img,2),size(img,3),5*n_te]);
    end
    
    
slc_set=1;    
    
for ii_slice=slc_set%1:24
    clear convGrappa_mag
    clear Img_SuperRes Img_lowRes
    
convGrappa_mag=img_ini(:,:,ii_slice,:);

complex_data_conv=reshape(convGrappa_mag,[size(convGrappa_mag,1),size(convGrappa_mag,2),size(convGrappa_mag,3),nrf_d,size(convGrappa_mag,4)/nrf_d]);
complex_data_conv=permute(complex_data_conv,[1 2 3 5 4]);
%%

PF = 6/8;   % for 960um PF is 7/8 partial fourier
ExtraCropPElines = 1; % number of extra k-line to not trust next to p.f. zero filled space, due to grappa smoothing around this region. 
POCS_flag = 0; % do POCS for p.f.
avgs = 1;
resolution =1.00;% 1.00;
%% complex_data: x,y,z,diffdir,RFenc,
complex_data = complex_data_conv;


% clear phase_data mag_data
% use phase of filtered low Res image as estimate for background phase and remove it. (reasonable for low bvalues)
Img_lowRes = permute(RealDiffusion_lowRes(permute(complex_data,[2,1,3,4,5]),PF,ExtraCropPElines,POCS_flag),[2,1,3,4,5]);

%%
if nrf_d==10
    [n1,n2,nz,ndir0,nrf]=size(Img_lowRes);
    Img_lowRes = reshape(permute(reshape(Img_lowRes,[n1,n2,nz,2,ndir0/2,nrf]),[1,2,3,5,6,4]),[n1,n2,nz,ndir0/2,nrf*2]);
end
[nx,ny,nz_l,ndir,nrf_d]=size(Img_lowRes);
%% load slice profile 


nrf=nrf_d;

mxy=MOut;
z=zOut;
B1=B1map_de;
clear B1map B1map_intp B1map phi_m
% SliderShift = [0 0 0 0 0];
SliderShift = zeros(1, nrf);
SidelobesSlices = 6; 
CropEdge = 1;


[A_encoding,psf_matrix] = B1corr_SliderPSFmatrix_pShift_weightDownEdges(mxy,z,nz_l,nrf,SliderShift,SidelobesSlices,CropEdge,B1,B1_index);
%%
for ii=1:n_te
A_encoding_b((1:5)+(ii-1)*5,(1:5)+(ii-1)*5,:,:)=A_encoding;
end
for ii=1:5
    phi_m((1:n_te)+(ii-1)*n_te,(1:v_num)+(ii-1)*v_num)=V(:,1:v_num);
end

index_x=[1:5:5*n_te,2:5:5*n_te,3:5:5*n_te,4:5:5*n_te,5:5:5*n_te];
%index_x=[1:6:30,2:6:30,3:6:30,4:6:30,5:6:30,6:6:30];
A_encoding_c=A_encoding_b(:,index_x,:,:);
img_sel=floor(1:5*n_te);a=zeros(5,n_te);a(img_sel)=1

img_b=permute(Img_lowRes,[1,2,3,5,4]);
img_b=reshape(img_b,[200,200,1,1,5*n_te]);


A_encoding_d=zeros(size(A_encoding_c));
img_c=zeros(size(img_b));
if 1 % 30x30
%A_encoding_d(img_sel,:,:,:)=A_encoding_c(img_sel,:,:,:);
A_encoding_d=A_encoding_c;
img_c(:,:,1,1,img_sel)=img_b(:,:,1,1,img_sel);
else
A_encoding_d=A_encoding_c(img_sel,:,:,:);
img_c=img_b(:,:,1,1,img_sel);
end


mask=zeros(size(A_encoding_d,2),1);
mask(img_sel)=1;
%mask(:)=1;
%%
clear A_encoding_e
for nn=1%:size(A_encoding_d,3)
    for mm=1%:size(A_encoding_d,4)
A_encoding_e(:,:,nn,mm)=A_encoding_d(:,:,nn,mm)*(phi_m);
    end
end
%%


%%
lambdaTikPercent=0.4;%0.5
Img_SuperRes = gSliderReconstruction_ite_tik_test11(img_c,A_encoding_e,lambdaTikPercent,mask);


PhaseFac=pi/2;

Img_SuperRes = real(Img_SuperRes*exp(1i*PhaseFac));
%%
Img_SuperRes = permute(Img_SuperRes,[2,1,3,4]);
%%
disp('Done');

img_recon(:,:,ii_slice,:)=Img_SuperRes;

ii_slice
end

%save_nii(make_nii(Img_SuperRes, [1.00 1.00 1.00], [], 16), ['DWIs.nii'])

%% mrf
img_save=img_recon;

img_recon=reshape(img_recon,[size(img_recon,1),size(img_recon,2),size(img_recon,3),v_num,5]);
img_recon=permute(img_recon,[1,2,5,3,4]);
img_recon=reshape(img_recon,[size(img_recon,1),size(img_recon,2),size(img_recon,4)*5,v_num]);
img_coef=img_recon;

for ii_slice=slc_set*5-4:slc_set*5%1:size(img_recon,3)
        
        n_v=2;
        temp_a=img_recon(:,:,ii_slice,:);
        temp_a=permute(reshape(temp_a,[size(img_recon,1)*size(img_recon,2),size(img_recon,4)]),[2,1]);
        temp_a=V(:,1:n_v)*temp_a(1:n_v,:);
        temp_a=permute(temp_a,[2,1]);
        temp_a=reshape(temp_a,[size(img_recon,1),size(img_recon,2),n_te]);
        img_rec(:,:,ii_slice,:)=temp_a;
        

ii_slice
end

%%%%%%%%%%%%%%%%
mask_air=1;
bg=1;
ed=n_te;
n_window=1;
ed=ed-n_window+1;
TR_num=ed-bg+1;
TE_sel=[1:5];

sli_sel=[slc_set*5-4:slc_set*5];
for ii_sli=1:length(sli_sel)

%% load
clear IMG I_test I_test_p

img_recon=img_rec;
%img_recon=reshape(img_recon,[210,220,5,5]);
%img_recon=permute(img_recon,[1,2,4,3]);
%img_recon(:,:,:,4:5)=0;
img_all=img_recon+eps;

img=squeeze(img_all(:,:,sli_sel(ii_sli),:)); % transverse
%img=squeeze(img_all(:,sli_sel(ii_sli),:,:)); % sag

clear I I_dic
%load ./dic/dic_test6_6te_10shots_tr3000
%load ./dic/dic_test5_10shots%T2SIM_TE20
load (dic_path)
I=I(:,:,:,:);
T1=0;

mask_img=ones(size(img,1),size(img,2));


%% delete phase
img=abs(img);
I=abs(I);




%% window with phase
for ii=1:ed
    %I_phase(:,:,:,ii)=sum((I(:,:,:,ii:ii+n_window-1)),4);
    %I(:,:,:,ii)=sum(abs(I(:,:,:,ii:ii+n_window-1)),4);
  
    
     img(:,:,ii)=sum(img(:,:,ii:ii+n_window-1),3);
     I(:,:,:,ii)=sum((I(:,:,:,ii:ii+n_window-1)),4);
end
I=(I(:,:,:,TE_sel));
I_test=(img(:,:,TE_sel));

I_test_p=I_test; % retain the proton density
image_size=size(img);% TR=10+4*rand(1,TR_num);



%%

T1_find=zeros(image_size(1),image_size(2));
T2_find=zeros(image_size(1),image_size(2));
p=zeros(image_size(1),image_size(2));
for m=1:image_size(1)
    for n=1:image_size(2)
I_test(m,n,:)=I_test_p(m,n,:)/(sum(abs(I_test_p(m,n,:)).^2))^0.5;
%I_test(m,n,:)=I_test_p(m,n,:)/(max(abs(I_test_p(m,n,:)).^2))^0.5;
    end
end
%for I_T1=1:size(T1,2)
for I_T1=1:size(T1,2)
for I_T2=1:size(T2,2)
        for I_OR=1:size(OR,2)
        I_dic(I_T1,I_T2,I_OR,:)=I(I_T1,I_T2,I_OR,:)/(sum(abs(I(I_T1,I_T2,I_OR,:)).^2))^0.5;
 %      I_dic(I_T1,I_T2,I_OR,:)=I(I_T1,I_T2,I_OR,:)/(max(abs(I(I_T1,I_T2,I_OR,:)).^2))^0.5;
        end
    end
end

% I_dic=I/max(I(:));
% I_test=I_test/max(I(:));

[T1_find,T2_find,OR_find,p,diff] = mrf_recog_max_dotproduct_img_su_or_diff_sm (I_test,I_test_p,I_dic,I,T1,T2,OR,mask_img);
p=p.*(abs(diff));
%p=p/max(abs(p(:)));
%%

%%
if 0
figure;colormap hot;

% subplot(2,2,1),imagesc(image(:,:,1));
% title('fully-sampled img')
subplot(2,2,1),imagesc(T1_find);
title('recognized T1 map')
colorbar ;caxis ([0,4000]);
subplot(2,2,2),imagesc(T2_find);
 colorbar ;caxis ([0,500]);
title('recognized T2 map')
subplot(2,2,3),imagesc(abs(p));
colorbar ;caxis ([0,1]);
title('recognized proton density map')
subplot(2,2,4),imagesc(OR_find);
colorbar ;
title('off-resonance map')
% figure,imshow(OR_find,[]);
end
%%


%plot_dic_test_or((img),(I_test),(I_dic),T1_find,T2_find,OR_find,T1,T2,OR);

%% mask air
if 1
mask_air=ones(image_size(1),image_size(2));
air_loc=find(p(:)<max(p(:))*5e-2);
mask_air(air_loc)=0;
mask_air=medfilt2(mask_air,[5,5]);
T2_find=T2_find.*mask_air;
mask_air_all(:,:,ii_sli)=mask_air;
% figure;colormap hot;
% 
% % subplot(2,2,1),imagesc(image(:,:,1));
% % title('fully-sampled img')
% subplot(2,2,1),imagesc(T1_find.*mask_air);
% title('recognized T1 map')
% colorbar ;caxis ([0,4000]);
% subplot(2,2,2),imagesc(T2_find.*mask_air);
%  colorbar ;caxis ([0,4000]);
% title('recognized T2 map')
% subplot(2,2,3),imagesc(abs(p.*mask_air));
% colorbar ;caxis ([0,1]);
% title('recognized proton density map')
% subplot(2,2,4),imagesc(OR_find);
% colorbar ;
% title('off-resonance map')
end

%%
if 0
    
xInput=150;
yInput=110;
T1_num=find(T1==T1_find(xInput,yInput));
    T2_num=find(T2==T2_find(xInput,yInput));
    OR_num=find(OR==OR_find(xInput,yInput));
    figure,plot(abs(squeeze(I_test(xInput,yInput,:))))
    set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
    figure
plot(squeeze(abs(I_test(xInput,yInput,:))))
        hold,plot(abs(squeeze(I_dic(T1_num,T2_num,OR_num,:))),'r','linewidth',4)
         set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
       % title(strcat('location:(',num2str(xInput),',',num2str(yInput),')T1:',num2str(T1_find(xInput,yInput)),'  T2:',num2str(T2_find(xInput,yInput)),'  OR:',num2str(OR_find(xInput,yInput)),'  red:recognized curve'))
end

%%

T2_map_all(:,:,ii_sli)=T2_find;
p_all(:,:,ii_sli)=p;
ii_sli
end 

T2_map_all=rot_mat(T2_map_all,180);
figure;colormap hot;
imagesc(cat(2,T2_map_all(:,:,1),T2_map_all(:,:,5)))
%imagesc(rot_mat((T2_map_all(:,:,1),3),180));
 colorbar ;caxis ([0,300]);
 axis image
 
 