function [img] = qslider_T2W_sim (T2_map_all,p_all,TE_set)


TE=TE_set;
FA=90*ones(1,length(TE));
FA_i=180*ones(1,length(TE));
TR=3500*ones(1,length(TE));
T1_set=1000;

for ii_rf=1:size(T2_map_all,4)
 for ii_z=1:size(T2_map_all,3)
for ii_x=1:size(T2_map_all,1)
    for ii_y=1:size(T2_map_all,2)     
            
            [F0_vector_out,Xi_F_out,Xi_Z_out]=ssfp_epg_cxz_20150916_se_rfp_ir (FA,FA_i,TR,TE,T1_set,T2_map_all(ii_x,ii_y,ii_z,ii_rf),0);
            img(ii_x,ii_y,ii_z,ii_rf,:)=(F0_vector_out(:));
        end
    end
ii_z
ii_rf
end
end
img=img.*repmat(p_all,[1,1,1,1,size(img,5)]);