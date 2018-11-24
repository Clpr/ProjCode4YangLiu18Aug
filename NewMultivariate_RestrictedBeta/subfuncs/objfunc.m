function [LogL,G_expand] = objfunc(X,L,I,T,K,Y,Z)
%{
objective function of GARCH-MIDAS (multivariate, restricted 2 par Beta weighting)

Input:
1. X [col-vec,len=3L+4]: in order: [Mu,Alpha,Beta,M,Theta1,...,ThetaL,W1_1,...,W1_L,W2_1,...,W2_L]
1. L [int]: number of macro-economic variables
1. I [int]: number of observations in each period (t)
1. T [int]: number of periods TO ESTIMATE (also freq of macro-economic variables) (the very first K lag obs cut)
1. K [int]: number of lag for each macro-variable
1. Y [mat,size=T*I]: observed 1-dim data to depart (high freq), time from 1 to T (old to new) and 1 to I (old to new)
1. Z [mat,size=(T+K)*L]: macro-economic vars (low-freq); note, additional K lags required; time from 1 to T (old to new)   
1. flagWeight [int]: flag for which kind of weight strategy to use
    1 for Beta weighting (unstricted, 1 parameter)
    2 for restricted (2 parameters)
    3 for ALMON (2 parameters)

Output:
1. LogL [col-vec,len=T*I]: log-likelihood of each r_{it}

%}
% 
% % test
% L = 5; I = 20; T = 30; K = 10;
% X = rand(3*L+4,1);
% Y = randn(T,I); Z = randn(T+K,L);


% Cap 0: Data allocation
	Mu = X(1); 
	Alpha = X(2); Beta = X(3);
	M = X(4);
	Theta = X(5:4+L); W1 = X((5+L):(4+2*L)); W2 = X((2*L+5):end);

	% validation for GARCH weighting (non-negative weights)
	if 1-Alpha-Beta<0
		LogL = -Inf;
		return;
	end
	


% Cap1: Weighting Strategy (2 par Restricted Beta weighting strategy)
	Weights = zeros(T*K,L);  % lagged-periods (t-1,...,t-K) for t=1,...,T * for each macro-var

	for idxT = T:-1:1
		% 1.1 power index
		seq = transpose(1:K);
		% 1.2 compute weight
		chk = (idxT-1)*K+1 : idxT*K ;  % short writing for indexing
		Weights(chk,:) = repmat((1-seq./K+10*eps),1,L);
        tmp_part2 = repmat(seq./K,1,L);
		Weights(chk,:) = Weights(chk,:) .^( repmat(W1',K,1) - 1 )  .*  tmp_part2.^( repmat(W2',K,1) - 1 )  ;
		Weights(chk,:) = Weights(chk,:) ./ repmat( sum(Weights(chk,:)), K,1) ;
	end
      
% Cap2: Estimate \tau_{t} (size=T*1) (long-term volatility)
	Tau = zeros(T,1);
	for idxT = T:-1:1
		chk_w = (idxT-1)*K+1 : idxT*K ;  % indexes to select weights from Weights
		chk_z = idxT : idxT+K-1 ; % indexes to select macro-vars from Z
		Tau(idxT) = exp( M + sum(Weights(chk_w,:).*Z(chk_z,:),1) * Theta );
	end

% Cap3: Estimate g_{it} (short-term volatility) (size=T*I,1)
	% an initial value of 1 for g_{1,1} is workable; it only affects the scaling of \beta estimate which we do not care
	G_expand = ones(T*I,1);
	Tau_expand = reshape( transpose(Tau*ones(1,I)) , I*T,1); % expand \tau to each i, for convenience of vetorized expr
	Y_expand = reshape(Y',T*I,1);
	for idxIT = 2:T*I
		G_expand(idxIT) = (1-Alpha-Beta) + Alpha*(Y_expand(idxIT-1) - Mu).^ 2 ./ Tau_expand(idxIT) + Beta*G_expand(idxIT-1);
	end
	% reshape G to matrix
% 	G = reshape(G_expand,T,I);


% Cap4: Compute quasi Log-likelihood
	Variance = abs(Tau_expand.*G_expand) ;
	if any(Variance)<0
		LogL = -Inf;
		return;
	end
	
	LogL = -0.5 * ( log(2*pi.*Variance) + (Y_expand-Mu).^2 ./ Variance );
    LogL = sum(LogL);
	return;
	
end  % function ends

