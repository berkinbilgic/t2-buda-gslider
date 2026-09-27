function Img_SuperRes =gSliderReconstruction_ite_tik_test10(Img_lowRes,A_recon,lambdaTikPercent,mask)
num_iter=10;
fista=1;
wavelet_flag=0;  
tik_flag=1;
mask_flag=1;

mask_coef_flag=1;
th_flag=1;

th_threshold=0.1e-3; %% 0.1e-3
wavelet_lambda=0.0008; % 0.0012
tol=0.01; %default: 0.01
tik=1; %0.4




[nx,ny,nz_l,ndir,nrf]=size(Img_lowRes);
n_basis = size(A_recon,1)/size(A_recon,2);
Img_SuperRes = zeros(nx,ny,nz_l*nrf/n_basis,ndir);
A_recon=A_recon(:,:,1,1);

      %%  tiknov reg
      if tik_flag
        A_recon = cat(1,A_recon,tik*eye(size(A_recon,2),size(A_recon,2)));
         mask((size(Img_lowRes,5)+1):(size(Img_lowRes,5)+size(A_recon,2)))=0;
        Img_lowRes(:,:,:,:,(size(Img_lowRes,5)+1):(size(Img_lowRes,5)+size(A_recon,2)))=0;
      end
      
%%        
        A_tik = InverseAmatrix(A_recon,lambdaTikPercent); %c
      % A_tik=pinv(A_recon);
%         vox=squeeze(Img_lowRes(x,y,:,:,:));
%         rhs=permute(vox, [3,1,2]);
%         rhs=reshape(rhs, size(rhs,1)*size(rhs,2),size(rhs,3));
%         res=A_tik *rhs;
%         Img_SuperRes(x,y,:,:)=res;
        for diff_dir = 1:ndir  
            img_b_ini= squeeze(Img_lowRes(:,:,1,diff_dir,:)); % 3rd: slice index
            img_size=size(img_b_ini); % 3rd: RF-TE dimension
            
            
            img_b_row=reshape(img_b_ini,[img_size(1)*img_size(2),img_size(3)]);
            
            img_b_row = permute(img_b_row, [2,1]); % rhs
            
            
            t_k=1;
          
            y_k= A_tik * img_b_row; 
            x_kneg1=y_k;
            W=Wavelet('Daubechies',4,4);
            wavWeight=0.4;
            %% iteration
for ii_iter=1:num_iter
    if fista
        t_kplus1 = (1 + sqrt(1 + 4 * t_k^2)) / 2;
        coef_kneg1 = -(t_k - 1) / t_kplus1; 
        coef_k = (t_kplus1 + t_k - 1) / t_kplus1; 
    else
        coef_k = 1;
        coef_kneg1 = 0;
        t_kplus1 = 1;
    end 
            img_coef_row=y_k;  
     %      img_coef_row=img_coef_row+A_tik*(img_b_row-A_recon*img_coef_row);
           img_coef_row=img_coef_row+A_tik*(img_b_row-A_recon*img_coef_row.*repmat(mask,[1,size(img_coef_row,2)]));
 %          img_coef_row=img_coef_row+A_tik*(img_b_row-A_recon*img_coef_row).*repmat(mask,[1,size(img_coef_row,2)]);
     
    
  %% mask
    if mask_flag
        img_coef_row=A_recon*img_coef_row;
        img_coef_row=img_coef_row.*repmat(mask,[1,size(img_coef_row,2)]);
        img_coef_row=A_tik*img_coef_row;       
    end
    
  %% mask_coef
    if mask_coef_flag
        n_te=size(img_coef_row,1)/5;
        n_bg=3;
        mask_coef=cat(2,[n_bg:n_te],[n_bg:n_te]+n_te,[n_bg:n_te]+2*n_te...
            ,[n_bg:n_te]+3*n_te,[n_bg:n_te]+4*n_te);
        img_coef_row(mask_coef,:)=0;
    end
             if th_flag
                 %  th_threshold=1e-2;
            decay=0.99;
            [U,S,V]=svd((img_coef_row),'econ');
            S=S-th_threshold*max(S(:));
            S(S(:)<0)=0;
            img_coef_row=(U*S*V');
            th_threshold=th_threshold*decay;
            end
         
            img_coef=permute(reshape(img_coef_row,[size(img_coef_row,1),img_size(1),img_size(2)]),[2,3,1]);
  %% wavelet   
    if  wavelet_flag  
            img_coef=zpad(img_coef,[512,512,size(img_coef,3)]); % to diadic
            img_coef=W*img_coef;
            img_coef=softThresh(img_coef,wavelet_lambda*max(abs(img_coef(:))));%0.002
           % img_coef=softThresh(img_coef,wavWeight);
            img_coef=W'*img_coef;
            img_coef=crop(img_coef,[img_size(1),img_size(2),size(img_coef_row,1)]);
    end      
  %% 
            if 0
            [U,S,V]=svd(diag(img_coef),'econ');
            S=S-threshold*max(S(:));
            S(S(:)<0)=0;
            img_coef=diag(U*S*V');
           % threshold=threshold*decay;
            end
            
             x_k=permute(reshape(img_coef,[img_size(1)*img_size(2),size(img_coef_row,1)]),[2,1]);
          %  x_k=img_coef_row;
            
          %% update
            update=rmse(x_k,x_kneg1);
            %update=rmse(x_k(1:2:45,:),A_tik(1:2:45,:)*img_b_row);
          %%  
            y_kplus1 = coef_k*x_k+coef_kneg1*x_kneg1;
            
            t_k=t_kplus1;
            y_k=y_kplus1;
            x_kneg1=x_k;
            if update<tol
                break
            end
 disp(['iteration: ', num2str(ii_iter), '  update: ', num2str(update), ' %   coef_k: ', num2str(coef_k), '   coef_k-1: ', num2str(coef_kneg1)])
end
            
            
            Img_SuperRes(:,:,:,diff_dir) = permute(reshape(x_k,[size(x_k,1),img_size(1),img_size(2)]),[2,3,1]);
        end       
    end
 %   disp(['SuperRes: nx=', num2str(x), ' / ', num2str(nx)])
 


function [rmse]=rmse(in,true)
rmse=100*norm(in(:)-true(:))/norm(true(:));
end

function x = softThresh(y,t)
% apply joint sparsity soft-thresholding 
absy = sqrt(sum(abs(y).^2,3));
unity = y./(repmat(absy,[1,1,size(y,3)])+eps);

res = absy-t;
res = (res + abs(res))/2;
x = unity.*repmat(res,[1,1,size(y,3)]);
end