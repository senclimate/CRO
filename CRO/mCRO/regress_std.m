function [b, bint, r] = regress_std(y, X, alpha)
    % REGRESS_STD performs linear regression with standardized variables,
    % automatically detecting and handling intercept terms.

    if nargin < 3
        alpha = 0.05;
    end

    % Copy original X for later use
    X_orig = X;

    % Detect intercept column (column with all ones)
    intercept_idx = find(all(X == 1));
    has_intercept = ~isempty(intercept_idx);

    % Remove intercept column for standardization
    X_noint = X;
    X_noint(:, intercept_idx) = [];  % if no intercept, does nothing

    % Standardize X (z-score), but skip std=0 cases
    mu_X = mean(X_noint, 1);
    sigma_X = std(X_noint, 0, 1);
    sigma_X(sigma_X == 0) = 1;  % avoid division by zero
    X_std = (X_noint - mu_X) ./ sigma_X;

    % Standardize y
    mu_y = mean(y);
    sigma_y = std(y);
    if sigma_y == 0
        sigma_y = 1;
    end
    y_std = (y - mu_y) / sigma_y;

    % Add intercept column back if it existed
    if has_intercept
        X_reg = [ones(size(X,1),1), X_std];
    else
        X_reg = X_std;
    end

    % Regress standardized data
    [b_std, bint_std, r_std] = regress(y_std, X_reg, alpha);

    % Recover coefficients in original scale
    b = zeros(size(X,2),1);
    bint = zeros(size(X,2),2);

    if has_intercept
        % Non-intercept coefficients
        non_idx = setdiff(1:size(X,2), intercept_idx);
        b(non_idx) = b_std(2:end) .* (sigma_y ./ sigma_X(:));
        bint(non_idx,:) = bint_std(2:end,:) .* (sigma_y ./ sigma_X(:));

        % Intercept
        b(intercept_idx) = mu_y - sum((mu_X ./ sigma_X) .* b_std(2:end)) * sigma_y;
        bint(intercept_idx,1) = mu_y - sum((mu_X ./ sigma_X) .* bint_std(2:end,2)) * sigma_y;
        bint(intercept_idx,2) = mu_y - sum((mu_X ./ sigma_X) .* bint_std(2:end,1)) * sigma_y;
    else
        b = b_std .* (sigma_y ./ sigma_X(:));
        bint = bint_std .* (sigma_y ./ sigma_X(:));
    end

    % Residuals
    y_hat = X_orig * b;
    r = y - y_hat;
end


% function [par, par_int, residual] = regress_std(Y, X)
% % Standarize
% mu_X=mean(X,1);
% sigma_X=std(X,0,1);
% 
% mu_Y=mean(Y, 1);
% sigma_Y=std(Y, 0, 1);
% 
% % Standardize
% X_std=(X-mu_X)./sigma_X;
% Y_std=(Y-mu_Y)./sigma_Y;
% 
% % Regression
% par=regress(Y_std,X_std);
% 
% % Recover to Origianl Parameter
% par=(sigma_Y./sigma_X).*par;
% 
% Y_fit = X * par;
% residual = Y - Y_fit;
% 
% end



