function H = win_cos(x1,x2,param)
% cosine window

    omega = pi/(param.domain(2)-param.domain(1));
    h1 = sin(omega*(x1(:)-param.domain(1)));
    h2 = sin(omega*(x2(:)-param.domain(1)));
    H = h1*h2';

end