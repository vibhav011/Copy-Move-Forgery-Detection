function T = ransac3(img, matches) 
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
        T = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 
        
        % plotting only - comment if not needed
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
        plot(xFit2, yFit2, 'r-', 'LineWidth', 2); % Plot fitted line.      
        
        hold off;
end
