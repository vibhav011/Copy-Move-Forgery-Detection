function [theta, alpha, p, q] = ransac4(matches, T1)
    k = size(matches,1);
    x1 = matches(:, 1, 2);  
    y1 = matches(:, 1, 1);
    x2 = matches(:, 2, 2);  
    y2 = matches(:, 2, 1);
    coefficients1 = polyfit(x1, y1, 1);
    coefficients2 = polyfit(x2, y2, 1);
    m1 = coefficients1(1,1);
    m2 = coefficients2(1,1);
    m = (m1-m2)/(1+m1*m2);
    theta = atan(m);
%         T = [cos(theta) -sin(theta); sin(theta) cos(theta)]';

    ex=(T1*[1;0]);
    ey=(T1*[0;1]);
    alpha=(norm(ex)+norm(ey))/2;
%         T = alpha*T;

    A = zeros(k, 2);
    B = zeros(k, 2);

    A(:,:) = matches(:, 1, :);
    B(:,:) = matches(:, 2, :);

    muA = mean(A, 1); 
    muB = mean(B, 1);

    p = muA';
    q = muB';

%     x0 = muA' - T * muB';

    if alpha > 1
        theta = -theta;
        alpha = 1/alpha;
        tem = q;
        q = p;
        p = tem;
    end
end
