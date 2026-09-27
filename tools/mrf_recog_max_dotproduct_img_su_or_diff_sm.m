function [T1_find,T2_find,OR_find,p,diff] = mrf_recog_max_dotproduct_img_su_or_diff_sm (I_test,I_test_p,I_dic,I,T1,T2,OR,mask_img)

image_size=size(I_test);

down_rate=1;
down_select=image_size(1)*image_size(2)/down_rate;



T1_find=zeros(image_size(1),image_size(2));
T2_find=zeros(image_size(1),image_size(2));
OR_find=zeros(image_size(1),image_size(2));
% p=zeros(image_size(1),image_size(2));
TR_num=size(I_test,3);
dic_find=zeros(image_size(1),image_size(2));

I_2d=reshape(I,[size(I_dic,1)*size(I_dic,2)*size(I_dic,3),size(I_dic,4)]);
I_test_p_2d=reshape(I_test_p,[size(I_test,1)*size(I_test,2),size(I_test,3)]);

I_dic_2d=reshape(I_dic,[size(I_dic,1)*size(I_dic,2)*size(I_dic,3),size(I_dic,4)]);
%%
%I_dic_2d=I_dic_2d';
I_dic_2d=I_dic_2d.';
%%
I_test_2d=reshape(I_test,[size(I_test,1)*size(I_test,2),size(I_test,3)]);



% clear I_test I_dic 
for down_i=1:down_rate

bg_p=(down_i-1)*down_select+1;
ed_p=down_i*down_select;
    
I_test_2d_temp=I_test_2d(bg_p:ed_p,:);
dic_test=I_test_2d_temp*I_dic_2d;


% matlabpool open 6
for ii=1:size(dic_test,1)
%     dic_find(ii)=find(dic_test(ii,:)==max(dic_test(ii,:)));
temp=find(dic_test(ii,:)==max(dic_test(ii,:)));

dic_find(ii + (down_i-1)*down_select )=temp(1);
diff(ii + (down_i-1)*down_select )=max(dic_test(ii ,:));

end

clear dic_test I_test_2d_temp
down_i
end

T1_find=T1(mod(mod(dic_find-1,size(T1,2)*size(T2,2))+1-1,size(T1,2))+1).*mask_img;
%T2_find=T2(floor((mod(dic_find-1,size(T1,2)*size(T2,2))+1-1)/size(T1,2))+1).*mask_img;
T2_find=T2(floor((mod(dic_find-1,size(T1,2)*size(T2,2))+1-1)/size(T1,2))+1).*mask_img;
OR_find=OR(floor((dic_find-1)/(size(T1,2)*size(T2,2)))+1).*mask_img;

I_dic_p_2d=I_2d(dic_find,:);  % pre
%I_dic_pp_2d=I_dic_2d';
%I_dic_p_2d=I_dic_pp_2d(dic_find,:);

%% calculate pd 
%p=abs((sum(I_test_p_2d.^2,2)).^0.5./((sum(I_dic_p_2d.^2,2)).^0.5)); % pre
p=(sum(abs(I_test_p_2d).^2,2)).^0.5./((sum(abs(I_dic_p_2d).^2,2)).^0.5);
%p=(sum(abs(I_test_p_2d),2))./((sum(abs(I_dic_p_2d),2)));

p=reshape(p,[image_size(1),image_size(2)]).*mask_img;
%%
diff=reshape(diff,[image_size(1),image_size(2)]);
% a=1; %for fun

     %  figure,imshow(abs(diff),[])
   
end
