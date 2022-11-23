function fisher = cal_CRB_of_RIS_cali_correlate_chan(permtrSet,risGeniePhsSet,Hcas)
% ***************************************
%  calculate FIM
%  author - Wei Zhang
%  input-    
%            permtrSet: permutation matrice
%            risGeniePhsSet: genie RIS phase shift 
%            Hcas: effective channel
%  output: 
%             fisher: FIM 
%copyright - CSRL@Fudan,2022/11/23
%  ************************************
Mr = size(Hcas,1);
[O,Mris,Nb,~] = size(permtrSet);
mrisNb = Mris*Nb;
mrNb = Mr*Nb;
mrNbO = mrNb*O;
tmp1 = risGeniePhsSet.';
g = tmp1(:);
tmp2 = kron(Hcas, eye(Nb));

tolNum1 = Mris*(Nb-1);
tolNum2 = Mr*Mris;
tolNum = tolNum1+2*tolNum2;

%% muDerPhi
muDerEta = zeros(mrNbO,tolNum);
muDerEta(:,1:tolNum1) = cal_muDerPhi();
%% muDerHbreve
muDerEta(:,tolNum1+(1:tolNum2)) = cal_muDerh();
muDerEta(:,end-tolNum2+1:end) = 1j*muDerEta(:,tolNum1+(1:tolNum2));
%% covariance
% Roff = kron(eye(Mr),ones(Nb));
% Rz = (kron(ones(O),Roff)+kron(eye(O),eye(mrNb)));
% fisher = muDerEta'/Rz*muDerEta;
fisher = cal_fisher();
fisher = real(fisher);
    function y = cal_muDerPhi()
        y = zeros(mrNbO,Mris*(Nb-1));
        Pi = zeros(mrisNb*O,mrisNb);
        for oo = 1 : O
            for mm = 1:Mris
                Pi(mrisNb*(oo-1)+(mm-1)*Nb+(1:Nb),(mm-1)*Nb+(1:Nb)) = squeeze(permtrSet(oo,mm,:,:));
            end
        end
        tmp2Tl = kron(eye(O),tmp2);        
        tmp3 = 1j*(tmp2Tl*Pi)*diag(g);
        for mm = 1 : Mris
            y(:,(mm-1)*(Nb-1)+(1:Nb-1)) = tmp3(:,(mm-1)*Nb+(2:Nb));
        end
    end
    function y = cal_muDerh()
        y = zeros(mrNbO,tolNum2);
        for oo = 1 : O
            for ii = 1 : Mris
                for jj = 1 : Mr
                    Pioj = squeeze(permtrSet(oo,ii,:,:));
                    gj = tmp1(:,ii);
                    y((oo-1)*mrNb+(jj-1)*Nb+(1:Nb),(jj-1)*Mris+ii) = Pioj*gj;
                end
            end
        end
    end
    function y = cal_fisher()
        tmpp1 = muDerEta';
        tmpp2 = 0;
        
        for oo = 1 : O
            tmpp2 = tmpp2 + tmpp1(:,mrNb*(oo-1)+(1:mrNb));
        end
        
        tmpp6 = 0;
        for mm = 1 : Mr
            tmpp3 = tmpp2(:,(mm-1)*Nb+(1:Nb));
            tmpp4 = sum(tmpp3,2);
            tmpp5 = tmpp4*tmpp4';
            tmpp6 = tmpp6 + tmpp5;
        end
        y = tmpp1*tmpp1'-1/(O*Nb+1)*tmpp6;
    end
end