% Code created by Loïc Marrec

function y = sigm(x, theta, n)

    y = 1./(1+(x./theta).^n);
    
end
