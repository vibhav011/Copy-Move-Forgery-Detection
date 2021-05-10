function [T, x0] = ransac3(img, matches, T1)
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
        T = [cos(theta) -sin(theta); sin(theta) cos(theta)]';
        
        ex=(T1*[1;0]);
        ey=(T1*[0;1]);
        alpha=(norm(ex)+norm(ey))/2;
        T = alpha*T;
        
        A = zeros(k, 2);
        B = zeros(k, 2);
        
        A(:,:) = matches(:, 1, :);
        B(:,:) = matches(:, 2, :);
        
        muA = mean(A, 1); 
        muB = mean(B, 1);
        
        x0 = muA' - T * muB';
        
        if alpha > 1
            T = T' ./ (alpha^2);
            x0 = -T*x0;
        end
        
%         figure(56);
%         imshow(img);
%         drawnow;
%         hold on;
%         scatter(x1, y1, 'r.');
%         hold on;
%         scatter(x2, y2, 'b.');
%         hold off;
%         
%         plotting only - comment if not needed
        figure(8);
        imshow(img);
        drawnow;
        hold on;
        
        xFit1 = linspace(min(x1), max(x1), 1000);
        yFit1 = polyval(coefficients1 , xFit1);
        plot(xFit1, yFit1, 'r-', 'LineWidth', 2); % Plot fitted line.
        
        hold on;
        
        xFit2 = linspace(min(x2), max(x2), 1000);
        yFit2 = polyval(coefficients2 , xFit2);
        plot(xFit2, yFit2, 'b-', 'LineWidth', 2); % Plot fitted line.      
        
        hold off;
end
