function [ Weights ] = GetWeights_ExpAlmon( W1,W2,T,K,L )
%{
Computes weights matrix (Normalized Exponential Almon Weighting Strategy, 2 par)

Input:
    1. W1 [mat,size=L*1]: first par vector
    1. W2 [mat,size=L*1]: second par vector
    1. T [int]: periods
    1. K [int]: lag
    1. L [int]: number of macro vars

Output:
    1. Weights [mat,size=TK*L]: weighting matrix, K,...,K (total T times) by row
%}

        Weights = zeros(1*K,L);  % lagged-periods (t-1,...,t-K) for t=1,...,T * for each macro-var
	for idxT = T:-1:1
		% 1.1 power index
		seq = transpose(1:K);
		% 1.2 compute weight
		chk = (idxT-1)*K+1 : idxT*K ;  % short writing for indexing
		Weights(chk,:) = exp( seq*W1' + (seq.^2)*W2' )  ;
		Weights(chk,:) = Weights(chk,:) ./ repmat( sum(Weights(chk,:)), K,1) ;
    end

end

