function [Fcol,column_ix] = bestcolumn_2(Y,Rmax)
% Select sub-set of columns close to the best subpsace of R-best column
% space of Y
% Phan Anh-Huy 2021

normcol = sqrt(sum(Y.^2));

Ynorm = Y * diag(1./normcol);

sz = size(Y);

% SVD
[u, ~, ~] = svds(Y, 1*Rmax);


% Find the first column
%Q = Ynorm' * (eye(sz(1)) - u * u') * Ynorm + speye(sz(2));
%dQ2 = diag(Q)';
dQ = sum((Ynorm - u * (u' * Ynorm)).^2,1); % +1

% 
column_ix = [];
Fcol = [];

for no_currselected = 1:Rmax
    fprintf('Select Component %d ', no_currselected)

    QK = dQ;
    if (no_currselected > 1)
        QK = dQ + sum((FcolU' * Ynorm).^2,1);
    end

    [~,ix] = min(QK);

    fprintf('Column index %d\n',ix)

    % Indices of Selected columns
    column_ix = [column_ix ix];

    % new selected column
    Colnew = Y(:,ix);
    
    % Concatenate columns
    Fcol = [Fcol Colnew];

    % find subspace of selected columns
    [FcolU, r] = qr(Fcol,0);

end