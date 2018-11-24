function [ Weights ] = GetWeights_RestrictBeta( W1,W2,T,K,L )
%{
Computes weights matrix (Restricted Beta Weighting Strategy, 2 par)

Input:
    1. W1 [mat,size=L*1]: first par vector
    1. W2 [mat,size=L*1]: second par vector
    1. T [int]: periods
    1. K [int]: lag
    1. L [int]: number of macro vars

Output:
    1. Weights [mat,size=TK*L]: weighting matrix, K,...,K (total T times) by row
%}

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



end