function [Phi,Phio,lambda,zeta, U,S,V,U_r,S_r,V_r,W_r,W,Atilde1,Atilde,r] = DMD(x,thresh)
    X11 = x(:, 1:end-1);
    X22 = x(:, 2:end);
    m = size(X11, 2);

%% DMD
[U, S, V] = svd(X11, 'econ');
 Atilde1 = U' * X22 * V / S;
    [W, D] = eig(Atilde1); 
    Phio = X22 * V / S * W;
 
    r = length(find(diag(S)>thresh));
   %truncation to rank r
    U_r = U(:, 1:r);
    S_r = S(1:r, 1:r);
    V_r = V(:, 1:r);
    Atilde = U_r' * X22 * V_r / S_r;   
    [W_r, D] = eig(Atilde);
    Phi = X22 * V_r / S_r * W_r;
zeta = pinv(Phi)*X11;

lambda = diag(D);

