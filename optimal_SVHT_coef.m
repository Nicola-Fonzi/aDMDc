function coef = optimal_SVHT_coef(beta,Sig)

% function omega = optimal_SVHT_coef(beta, sigma_known)
%
% Coefficient determining optimal location of Hard Threshold for Matrix
% Denoising by Singular Values Hard Thresholding when noise level is known or
% unknown.  
%
% See D. L. Donoho and M. Gavish, "The Optimal Hard Threshold for Singular
% Values is 4/sqrt(3)", http://arxiv.org/abs/1305.5870
%
% IN: 
%    beta: aspect ratio m/n of the matrix to be denoised, 0<beta<=1. 
%          beta may be a vector 
% 
% OUT: 
%    coef:   optimal location of hard threshold, up the median data singular
%            value (sigma unknown) or up to sigma*sqrt(n) (sigma known); 
%            a vector of the same dimension as beta, where coef(i) is the 
%            coefficient correcponding to beta(i)
%
% Usage in unknown noise level:
%
%   Given an m-by-n matrix Y known to be low rank and observed in white
%   noise with mean zero and unknown variance, form a denoised matrix 
%   Xhat by:
%  
%   [U D V] = svd(Y); 
%   y = diag(Y); 
%   y( y < (optimal_SVHT_coef_sigma_unknown(m/n,0) * median(y)) ) = 0; 
%   Xhat = U * diag(y) * V';

coef = optimal_SVHT_coef_sigma_unknown(beta);
coef = coef*median(diag(Sig));

end

function lambda_star = optimal_SVHT_coef_sigma_known(beta)
    assert(all(beta>0));
    assert(all(beta<=1));
    assert(prod(size(beta)) == length(beta)); % beta must be a vector
    
    w = (8 * beta) ./ (beta + 1 + sqrt(beta.^2 + 14 * beta +1)); 
    lambda_star = sqrt(2 * (beta + 1) + w);
end

function omega = optimal_SVHT_coef_sigma_unknown(beta)
    warning('off','MATLAB:quadl:MinStepSize')
    assert(all(beta>0));
    assert(all(beta<=1));
    assert(prod(size(beta)) == length(beta)); % beta must be a vector
    
    coef = optimal_SVHT_coef_sigma_known(beta);

    MPmedian = zeros(size(beta));
    for i=1:length(beta)
        MPmedian(i) = MedianMarcenkoPastur(beta(i));
    end

    omega = coef ./ sqrt(MPmedian);
end




function med = MedianMarcenkoPastur(beta)
    MarPas = @(x) 1-incMarPas(x,beta,0);
    lobnd = (1 - sqrt(beta))^2;
    hibnd = (1 + sqrt(beta))^2;
    change = 1;
    while change && (hibnd - lobnd > .001)
      change = 0;
      x = linspace(lobnd,hibnd,5);
      for i=1:length(x)
          y(i) = MarPas(x(i));
      end
      if any(y < 0.5)
         lobnd = max(x(y < 0.5));
         change = 1;
      end
      if any(y > 0.5)
         hibnd = min(x(y > 0.5));
         change = 1;
      end
    end
    med = (hibnd+lobnd)./2;
end

function I = incMarPas(x0,beta,gamma)
    if beta > 1
        error('betaBeyond');
    end
    topSpec = (1 + sqrt(beta))^2;
    botSpec = (1 - sqrt(beta))^2;
    MarPas = @(x) IfElse((topSpec-x).*(x-botSpec) >0, ...
                         sqrt((topSpec-x).*(x-botSpec))./(beta.* x)./(2 .* pi), ...
                         0);
    if gamma ~= 0
       fun = @(x) (x.^gamma .* MarPas(x));
    else
       fun = @(x) MarPas(x);
    end
    I = quadl(fun,x0,topSpec);
    
    function y=IfElse(Q,point,counterPoint)
        y = point;
        if any(~Q)
            if length(counterPoint) == 1
                counterPoint = ones(size(Q)).*counterPoint;
            end
            y(~Q) = counterPoint(~Q);
        end
        
    end
end

