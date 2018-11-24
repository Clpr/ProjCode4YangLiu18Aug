clear;
clc;
%% ============= data reading ================
tab_Y = readtable('data/FX_JPY.csv');
tab_Z = readtable('data/MacroDat.csv');
tab_Z.nvix = tab_Z.nvix / 100;
% join(tab_Y,tab_Z);  % (left) join tables by common cols, requires high-freq data input as the 1st par
%% ============= PARMS DEF ================
% number of macro-vars
L = 18; 
% number of obs in one period
I = 20; 
% max lag
K = 12;
% number of periods (low-freq) (partially depends on K)
T = 179; 
%% ============= DATA PROCESS ================
% high-freq data matrix (size=T*I)
Y = zeros(T,I);
% low-freq data (macro var) (size=(T+K)*L)
Z = zeros(T+K,L);
% ----------------------------
% randomly delete obs in Y to make the number of obs in one period to be uniform (20)
counter = 1; seq_year = sort(unique(tab_Y.Year));
for idxY = 1:length(seq_year)
   for idxM = 1:12
      % jump from the last period, the March 2016 data are reserved for forecast
     if seq_year(idxY)==2016 && idxM==3
         break;
     end
     % select data of the specific month in specific year
     tmp = tab_Y.Y( (tab_Y.Year==seq_year(idxY))&(tab_Y.Month==idxM) );
     % real number of obs in this month
     tmp_count = length(tmp);
     % how many obs to delete, and randomly select obs from tmp
     tmp_delnum = tmp_count - 20;
     if tmp_delnum > 0
        rng(idxY*idxM); % constant seeds for repeat-ability
        tmp_del = randsample(tmp_count,tmp_delnum);
        tmp(tmp_del) = [];  % delete elements
     end
     Y(counter,:) = tmp';
     % counter + 1
     counter = counter + 1;
   end
end
% clear temp vars
clear counter seq_year idxY idxM tmp tmp_count tmp_delnum tmp_del;
% ----------------------------
% fill Z matrix (macro-vars)
tab_Z2 = tab_Z; tab_Z2.Year = []; tab_Z2.Month = [];
Z = table2array( tab_Z2 );
% ----------------------------
% cut useless part of Y (high-freq) and Z, from the very beginning (because of the lags of Z)
Z = Z(end-T-K+1:end,:);
Y = Y(end-T+1:end,:);

% use March 2016 data (the last period) as forecast dataset
fore_real_Y = table2array(  tab_Y( (tab_Y.Year==2016)&(tab_Y.Month==3)  , 4 )  );  % size=I*1 (but I >= 20)
fore_real_Z = table2array(  tab_Z( (tab_Z.Year==2016)&(tab_Z.Month==3) ,3:end )  );


%% =============== Estimation ==================
% preprare initial guess of parameters (length=2L+4)
% X = [Mu,Alpha,Beta,M,Theta1,...,ThetaL,W1_1,...,W1_L,W2_1,...,W2_L]
X0 = [ mean(mean(Y)); 
       0.01; 
       0.85;
       0.1;
       0.3*ones(L,1); 
       3*ones(L,1);
       3*ones(L,1)  ];
% prepare parms for optimization toolkit
lb = [ -Inf;
       0;
       0;
       -Inf;
       -Inf(L,1);
       1.001*ones(L,1);
       1.001*ones(L,1); ];
ub = [ Inf;
       1;
       1;
       Inf;
       Inf(L,1);
       50*ones(L,1);
       50*ones(L,1)  ];
% an anonymous function for clear input to fmincon()
afunc = @(X) -1 * objfunc(X,L,I,T,K,Y,Z)  ;
% settings
options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','MaxFunctionEvaluations',20000);
% options = optimoptions('fmincon','Algorithm','sqp','Display','iter','MaxFunctionEvaluations',20000);
%% estimation
[X_est,~,exit_flag,~,~,Jacob,Hessi] = fmincon(afunc,X0,[],[],[],[],lb,ub,[],options);
% display exit_flag (>0 means converged)
disp(exit_flag);
%% compute SE and p
CovMat = inv(Hessi);
SE = diag(CovMat); SE(SE<0) = NaN; SE = sqrt(SE);
Zstat = X_est ./ SE;
pVal = 0.5 * erfc(0.7071 * abs (Zstat)) * 2;

% distributes estimated parameters
est_Mu = X_est(1);
est_Alpha = X_est(2);
est_Beta = X_est(3);
est_M = X_est(4);
est_Theta = X_est(5:4+L);
est_W1 = X_est((5+L):(2*L+4));
est_W2 = X_est((2*L+5):end);

% result demo
Namelist = cell(3*L+4,1);
Namelist(1:4) = {'Mu';'Alpha';'Beta';'M'};
Znamelist = tab_Z2.Properties.VariableNames;
for idx = 1:L
    Namelist{4+idx} = ['Theta_',Znamelist{idx}];
    Namelist{4+L+idx} = ['W1_',Znamelist{idx}];
    Namelist{4+2*L+idx} = ['W2_',Znamelist{idx}];
end
Significance = pVal<0.05;
tab_result = table(Namelist, X_est, SE, pVal,Significance  )
% output
writetable(tab_result,'output/result.csv');


%% =============== FORECAST ==================
% select/unpackage estimated pars
est_Mu = X_est(1); 
est_Alpha = X_est(2); 
est_Beta = X_est(3);
est_M = X_est(4);
est_Theta = X_est(5:4+L); 
est_W1 = X_est((5+L):(2*L+4));
est_W2 = X_est((2*L+5):end);

% forecast daily r in March 2016
fore_est_Y = zeros( length(fore_real_Y) , 1 );
fore_est_g = zeros( length(fore_real_Y) , 1 );
fore_est_tau = zeros( length(fore_real_Y) , 1 );
for idx = 1:length(fore_est_Y)
    % compute weights (Beta unrestricted)
    Weights = zeros(1*K,L);  % lagged-periods (t-1,...,t-K) for t=1,...,T * for each macro-var
	for idxT = 1:-1:1
		% 1.1 power index
		seq = transpose(1:K);
		% 1.2 compute weight
		chk = (idxT-1)*K+1 : idxT*K ;  % short writing for indexing
		Weights(chk,:) = repmat((1-seq./K+10*eps),1,L);
        tmp_part2 = repmat(seq./K,1,L);
		Weights(chk,:) = Weights(chk,:) .^( repmat(est_W1',K,1) - 1 )  .*  tmp_part2.^( repmat(est_W2',K,1) - 1 )  ;
		Weights(chk,:) = Weights(chk,:) ./ repmat( sum(Weights(chk,:)), K,1) ;
    end
    % forecast \tau & g
    fore_est_tau(idx) = exp( est_M + sum(Weights.*fore_real_Z,1) * est_Theta );
    if idx == 1
        [~,last_g] = objfunc(X_est,L,I,T,K,Y,Z);
        last_g = last_g(end);
        fore_est_g(idx) = (1-est_Alpha-est_Beta) + est_Alpha*(Y(end) - est_Mu )^2 / fore_est_tau(idx) + est_Beta*last_g;
    else
        fore_est_g(idx) = (1-est_Alpha-est_Beta) + est_Alpha*(fore_est_Y(idx-1) - est_Mu )^2 / fore_est_tau(idx) + est_Beta*fore_est_g(idx-1);
    end
    % forecase r (Y)
    fore_est_Y(idx) =  est_Mu + sqrt( fore_est_tau(idx) * fore_est_g(idx) );
end
% sum sqaured forecast error
SQE = sum( (fore_est_Y - fore_real_Y).^2 );


%% ================= IN-SAMPLE FORECAST OF \tau & g ==============
% NOTE: compute in-sample (e.g. 2013-2016) \tau & g series;
% it is an in-sample forecast, obtaining two series from 2013 to 2016 rather than a one-step-forward forecast (e.g. series in 2017)
% NOTE: 计算【估计时使用的数据，如2013-2016】的长期波动序列\tau和短期波动序列g；
% 注意，是样本内预测，得到的是2013-2016年的两条序列，而不是向前一步如2017年的两条波动序列
% 这部分的计算实际上就是objfunc()的算法，只不过没有计算最后的似然值而已；如果把数据和估完的参数输进objfunc()然后修改返回值让其返回计算过程中的tau和g，结果是完全一致的
% data used: Y, Z
% parameter used: est_X
    % 1. get estimated parameters --> 不同数目参数的模型记得要修改取用估计结果的语句，直接复制粘贴上面的就行
    est_Mu = X_est(1); 
    est_Alpha = X_est(2); 
    est_Beta = X_est(3);
    est_M = X_est(4);
    est_Theta = X_est(5:4+L); 
    est_W1 = X_est((5+L):(2*L+4));
    est_W2 = X_est((2*L+5):end);
    % 2. get weights （wrote a sub-function to keep compat)  --> 对于不同的权重方案，只需要修改底下调用的权重计算方法就行，其余再往下完全一致
    Weights = GetWeights_RestrictBeta( est_W1,est_W2,T,K,L );

    % 3. get long-term volatility \tau, 'est_Tau' is the output
        % NOTE: this section is (nearly) directly copied from the objfunc() function, you can refer to its definition to look into details about how i organize the algorithm & data structures
    est_Tau = zeros(T,1);
	for idxT = T:-1:1
		chk_w = (idxT-1)*K+1 : idxT*K ;  % indexes to select weights from Weights
		chk_z = idxT : idxT+K-1 ; % indexes to select macro-vars from Z
		est_Tau(idxT) = exp( est_M + sum(Weights(chk_w,:).*Z(chk_z,:),1) * est_Theta );
    end
    % 4. get short-term volatility g, 'est_g' is the output, equivalent to 'G_expand' in objfunc()
        % NOTE: we select g0 = 1 as the start point
    est_g = ones(T*I,1);
	Tau_expand = reshape( transpose(est_Tau*ones(1,I)) , I*T,1); % expand \tau to each i, for convenience of vetorized expr
	Y_expand = reshape(Y',T*I,1);
	for idxIT = 2:T*I
		est_g(idxIT) = (1-est_Alpha-est_Beta) + est_Alpha*(Y_expand(idxIT-1) - est_Mu).^ 2 ./ Tau_expand(idxIT) + est_Beta*est_g(idxIT-1);
    end

    % 5. finally, do plotting to demo the results
        % NOTE: cauz the two series have diff freq, use two panes; you can modify it as you like
    figure();
        % 5.1 low-freq volatility, \tau
        subplot(1,2,1);
        plot(est_Tau);title('In-sample $\tau_{t}$, the long-term volatility','Interpreter','LaTex');
        xlabel('Low-Freq Time: t'); ylabel('$\tau_{t}$','Interpreter','LaTex'); grid on; xlim([0,length(est_Tau)+1]);
        % 5.2 high-freq volatility, g
        subplot(1,2,2);
        plot(est_g);title('In-sample $g_{i,t}$, the short-term volatility','Interpreter','LaTex');
        xlabel('High-Freq Time: i'); ylabel('$g_{i,t}$','Interpreter','LaTex'); grid on; xlim([0,length(est_g)+1]);
    












