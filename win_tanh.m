function H = win_tanh(x1,x2,param,gamma)
% cosine window

    rho = 1/(param.domain(2)-param.domain(1));
    h1 = tanh_win(rho*(x1(:) - param.domain(1)),gamma);
    h2 = tanh_win(rho*(x2(:) - param.domain(1)),gamma);
    H = h1*h2';

end

function h = tanh_win(x,gamma)
    h = zeros(size(x));
    I = x <= .5;
    h(I) = tanh(gamma*x(I));
    h(~I) = tanh(gamma*(1-x(~I)));
end