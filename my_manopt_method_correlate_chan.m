function [fbreve,cost] = my_manopt_method_correlate_chan(hbreve,E,fbreve,Mr,Nb,O)
%% manifold optimization
mrNb = Mr*Nb;
oneVec = ones(Nb,1);
[~,len] = size(E);
constMtr2 = E'- cal_constMtr2();
iterNum = 2e2;
costFunc = zeros(iterNum,1);
[costFunc(1)] = cost_func(fbreve);
initT = 0;
for ii = 2 : iterNum
    %%
    if ii == 2
        oldGrad = cost_func_der(fbreve);
        dir = -oldGrad; % update direction
    end
    df0 = real(oldGrad'*dir);
    [stepSize,newf,initT] = back_tracking_line_search(fbreve,df0,dir,initT);
    costFunc(ii) = cost_func(newf);
    dirTransp = proj(newf,dir);
    %%
    newGrad = cost_func_der(newf);
    betaRcg = norm(newGrad)^2/norm(oldGrad)^2;
    dir = -newGrad+betaRcg*dirTransp;
    %%
    fbreve = newf;
    oldGrad = newGrad;
    if ii > 10 &&  (abs(costFunc(ii)-costFunc(ii-10))/costFunc(ii-10) < 1e-4)
        break;
    end
end
cost = costFunc(1:ii);
% figure; semilogy(cost)

    function y = cost_func(fbreve)
        yTmp = hbreve-E*fbreve;
%         tic;
        y = yTmp'*yTmp - block_mtr_mul(yTmp);
%         toc;
    end
    function y = cost_func_der(fbreve)
%         tic;
        ey = -constMtr2*(hbreve-E*fbreve);
%         toc;
        y = proj(fbreve,ey);
    end
    function y = retr(z, v, t)
        if nargin <= 2
            t = 1.0;
        end
        y = z+t*v;
        y = y ./ abs(y);
    end
    function y = proj(z,u)
        y = u - real( conj(u) .* z ) .* z;
    end
    function [stepSize,nextf,initT] = back_tracking_line_search(f,df0,direction,initT)
        alpha = 0.5;
        beta = 0.5;
        dirNorm = norm(direction);
        if initT == 0
            t = 1/dirNorm;
        else
            t = initT;
        end
        tmp2 = cost_func(f);
        costEva = 1;
        while 1
            nextf = retr(f, direction, t);
            tmp1 = cost_func(nextf);
            if tmp1 > tmp2 + alpha*t*df0
                if t < 1e-5
                    nextf = f;
                    break;
                end
                t = beta*t;
            else
                break;
            end
            costEva = costEva+1;
        end
        stepSize = t*dirNorm;
        switch costEva
            case 1
                % If things go well, push your luck.
                initT = 2*t;
            case 2
                % If things go smoothly, try to keep pace.
                initT = t;
            otherwise
                initT = 0;
        end
        
    end
    function y = block_mtr_mul(yTmp)
        
        tmp1 = zeros(mrNb,1);
        for oo = 1 : O
            tmp1 = tmp1 + yTmp((oo-1)*mrNb+(1:mrNb));
        end
        tmp2 = 0;
        for mm = 1 : Mr
            tmp3 = tmp1((mm-1)*Nb+(1:Nb));
            tmp4 = oneVec'*tmp3;
            tmp2 = tmp2 + tmp4'*tmp4;
        end
        y = 1/(1+O*Nb)*tmp2;
        
    end
    function y = cal_constMtr2()
        tmp1 = E';
        tmp2 = zeros(len,mrNb);
        for oo = 1 : O
            tmp2 = tmp2+tmp1(:,(oo-1)*mrNb+(1:mrNb));
        end
        
        for mm = 1 : Mr
            tmp3 = tmp2(:,(mm-1)*Nb+(1:Nb));
            tmp4 = (tmp3*oneVec)*oneVec';
            tmp2(:,(mm-1)*Nb+(1:Nb)) = tmp4;
        end
        
        y = 1/(1+O*Nb)*repmat(tmp2,1,O);
    end
end
