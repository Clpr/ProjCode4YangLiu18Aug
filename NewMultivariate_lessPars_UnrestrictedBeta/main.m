%% ============= data reading ================
clear;
clc;

tab_Y = readtable('data/FX_JPY.csv');
tab_Z = readtable('data/MacroDat.csv');
tab_Z.nvix = tab_Z.nvix / 100;
% join(tab_Y,tab_Z);  % (left) join tables by common cols, requires high-freq data input as the 1st par
% ============= PARMS DEF ================
% number of macro-vars
L = 18; 
% number of obs in one period
I = 20; 
% max lag
K = 12;
% number of periods (low-freq) (partially depends on K)
T = 180; 
% ============= DATA PROCESS ================
% high-freq data matrix (size=T*I)
Y = zeros(T,I);
% low-freq data (macro var) (size=(T+K)*L)
Z = zeros(T+K,L);
% ----------------------------
% randomly delete obs in Y to make the number of obs in one period to be uniform (20)
counter = 1; seq_year = sort(unique(tab_Y.Year));
for idxY = 1:length(seq_year)
   for idxM = 1:12
      % jump from the last period
     if seq_year(idxY)==2016 && idxM==4
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


%% =============== Estimation in complex area ==================
% % preprare initial guess of parameters (length=2L+4)
% % X = [Mu,Alpha,Beta,M,Theta1,...,ThetaL,W1,...,WL]
% X0 = [ mean(mean(Y)); 
%        0.05; 
%        0.5;
%        0.1;
%        0.1; 
%        3  ];
% % prepare parms for optimization toolkit
% lb = [ -Inf;
%        0;
%        0;
%        -Inf;
%        -Inf;
%        1.001 ];
% ub = [ Inf;
%        1;
%        1;
%        Inf;
%        Inf;
%        50 ];
% % an anonymous function for clear input to fmincon()
% afunc = @(X) call_obj(X,L,I,T,K,Y,Z)  ;
% % settings
% options = optimoptions('fmincon','GradObj','on','Algorithm','interior-point','Display','iter');
% % estimation
% [X_est,~,exit_flag,~,~,Jacob,Hessi] = fmincon(afunc,X0,[],[],[],[],lb,ub,[],options);
% % display exit_flag (>0 means converged)
% disp(exit_flag);
%% =========== ESTIMATION only in real area ===================
% preprare initial guess of parameters (length=2L+4)
% X = [Mu,Alpha,Beta,M,Theta1,...,ThetaL,W1,...,WL]
% X0 = [ mean(mean(Y)); 
       % 0.2; 
       % 0.5;
       % 0.1;
       % 0.3; 
       % 5  ];
   X0 = [ mean(mean(Y)); 
       0.5; 
       0.3;
       0.1;
       0.3; 
       2  ];
% prepare parms for optimization toolkit
lb = [ -Inf;
       0;
       0;
       -Inf;
       -Inf;
       1.001 ];
ub = [ Inf;
       1;
       1;
       Inf;
       Inf;
       50 ];
afunc2 = @(X) -1 * sum( objfunc(X,L,I,T,K,Y,Z) )  ;
options2 = optimoptions('fmincon','Algorithm','interior-point','Display','iter','MaxFunctionEvaluations',5000);
[X_est,~,exit_flag,~,~,Jacob,Hessi] = fmincon(afunc2,X0,[],[],[],[],lb,ub,[],options2);
disp(exit_flag);

%% compute SE and p
CovMat = inv(Hessi); 
SE = sqrt(diag(CovMat));
Zstat = X_est ./ SE;
pVal = 0.5 * erfc(0.7071 * abs (Zstat)) * 2;

% distributes estimated parameters
est_Mu = X_est(1);
est_Alpha = X_est(2);
est_Beta = X_est(3);
est_M = X_est(4);
est_Theta = X_est(5);
est_W = X_est(6);

% result demo
Namelist= {'Mu';'Alpha';'Beta';'M';'Theta';'W'};
tab_result = table(Namelist, X_est, SE, pVal )




%% griding at alpha and beta
% N = 100; 
% lin_Alpha = ones(N,1);
% lin_Beta = ones(N,1);
% n = 1;
% while n<=N
%     tmp_alpha = rand(); tmp_beta = rand();
%     if tmp_alpha + tmp_beta <0.9
%         lin_Alpha(n) = tmp_alpha;
%         lin_Beta(n) = tmp_beta;
%         n = n+1;
%     end
% end
% 
% grid_X = cell(N,1); grid_SE = cell(N,1); grid_pVal = cell(N,1);
% warning('off');
% for idxA = 1:N
%     % X = [Mu,Alpha,Beta,M,Theta1,...,ThetaL,W1,...,WL]
%     X0 = [ mean(mean(Y)); 
%    lin_Alpha(idxA); 
%    lin_Beta(idxA);
%    0.1;
%    0.3; 
%    5  ];
%     % prepare parms for optimization toolkit
%     lb = [ -Inf;
%            0;
%            0;
%            -Inf;
%            -Inf;
%            1.001 ];
%     ub = [ Inf;
%            1;
%            1;
%            Inf;
%            Inf;
%            50 ];
%     afunc2 = @(X) -1 * sum( objfunc(X,L,I,T,K,Y,Z) ) / (1E+3) ;
%     options2 = optimoptions('fmincon','GradObj','off','Algorithm','interior-point','Display','off');
%     [X_est,~,exit_flag,~,~,Jacob,Hessi] = fmincon(afunc2,X0,[],[],[],[],lb,ub,[],options2);
%     CovMat = inv(Hessi); CovMat(CovMat<0)=NaN;
%     SE = sqrt(diag(CovMat));
%     Zstat = X_est ./ SE;
%     pVal = 0.5 * erfc(0.7071 * abs (Zstat)) * 2;
%     % save results
%     grid_X{idxA} = X_est;
%     grid_SE{idxA} = SE;
%     grid_pVal{idxA} = pVal;
% 
% 
% end
% warning('on');
% disp('Done.');
% 
% 
% 
% 
% tmp = table( lin_Alpha, lin_Beta, grid_SE, grid_pVal, grid_X );

















