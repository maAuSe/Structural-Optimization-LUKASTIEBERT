function [U,fK] = lsolver(K,P,save_factor,check_crpoint,fK)

%LSOLVER   Solves linear system of equations.
%
%   [U,fK] = lsolver(K,P,save_factor,check_crpoint,fK)
%   solves the linear system K*U=P.
%
%   K    			Coefficient matrix		[n * n]
%   P   			Right-hand side     	[n * 1]
%   save_factor		Save factorization		{true | false}
%   check_crpoint	Check critical points	{true | false}
%   fK  			Struct containing factorization of K		

% select linear solver
if nargin < 5 || isempty(fK) % factorization of K is not provided
    if and(~save_factor,~check_crpoint) % not interested in factorization or critical point check -> backslash
        fK.type = 'backslash';
    else % perform factorization
        if isequal(K,K') % symmetric system -> ldl
            fK.type = 'ldl';
            [fK.L,fK.D,fK.p]=ldl(K,'vector');
            if check_crpoint; fK.d = full(diag(fK.D)); end
        else
            fK.type = 'lu';
            [fK.L, fK.U, fK.p, fK.q, fK.R] = lu(K,'vector');
            if check_crpoint; fK.d = full(diag(fK.U)); end
        end
    end
end


% solve system
if ~isempty(P)
    switch fK.type
        case 'backslash'
            U = K\P;
        case 'ldl'
            U = zeros(size(P));
            U(fK.p,:) =  fK.L.'\(fK.D\(fK.L\P(fK.p,:)));
        case 'lu'
            U(fK.q,:) = fK.U\(fK.L\(fK.R(:,fK.p)\P));
    end
else
    U = [];
end

end
