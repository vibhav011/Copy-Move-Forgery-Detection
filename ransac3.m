function [T, x0] = ransac3(img, matches, T1)
        k = size(matches,1);
        % coordinates of each patch 
        x1 = matches(:, 1, 2);  
        y1 = matches(:, 1, 1);
        x2 = matches(:, 2, 2);  
        y2 = matches(:, 2, 1);
        
        % best fit line for the two patches 
        coefficients1 = polyfit(x1, y1, 1);
        coefficients2 = polyfit(x2, y2, 1);
        
        % slopes of the best fit lines 
        m1 = coefficients1(1,1);
        m2 = coefficients2(1,1);
        
        % calculating the angle between the best fit line of the patches 
        m = (m1-m2)/(1+m1*m2);
        theta = atan(m);
        
        % rotation matrix 
        T = [cos(theta) -sin(theta); sin(theta) cos(theta)]';
        
        % scaling factor decision 
        ex=(T1*[1;0]);
        ey=(T1*[0;1]);
        
        if abs(theta) < 0.0523
            alpha=(norm(ex)+norm(ey))/2;
        else
            alpha=1;
        end
        
        % scaled rotation matrix
        T = alpha*T;
        
        % imposing conversion from smaller patch to larger one in case of scaled images
        A = zeros(k, 2);
        B = zeros(k, 2);
        
        A(:,:) = matches(:, 1, :);
        B(:,:) = matches(:, 2, :);
        
        muA = mean(A, 1); 
        muB = mean(B, 1);
        
        x01 = muA' - T * muB';
        x02 = muB' - T * muA';
        
        tr1 = T*A' + x02;
        tr2 = T*B' + x01;
        
        x0 = x01;
        if sum(vecnorm(tr1-B')) < sum(vecnorm(tr2-A'))
            x0 = x02;
        end
        
        if alpha > 1 
            T = T' ./ (alpha^2) % if needed rotation matrix is inverted and scaled down instead of up
            x0 = -T*x0 % translation factor is made negative 
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
%         figure(8);
%         imshow(img);
%         drawnow;
%         hold on;
%         
%         xFit1 = linspace(min(x1), max(x1), 1000);
%         yFit1 = polyval(coefficients1 , xFit1);
%         plot(xFit1, yFit1, 'r-', 'LineWidth', 2); % Plot fitted line.
%         
%         hold on;
%         
%         xFit2 = linspace(min(x2), max(x2), 1000);
%         yFit2 = polyval(coefficients2 , xFit2);
%         plot(xFit2, yFit2, 'b-', 'LineWidth', 2); % Plot fitted line.      
%         
%         hold off;
end
