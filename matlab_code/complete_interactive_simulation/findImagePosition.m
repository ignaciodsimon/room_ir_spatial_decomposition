function imagePosition = findImagePosition(sourcePosition, boundaryLimits)

    Ax = boundaryLimits(1);
    Ay = boundaryLimits(2);
    Bx = boundaryLimits(3);
    By = boundaryLimits(4);
    Px = sourcePosition(1);
    Py = sourcePosition(2);
% 
%     Ax = 0;
%     Ay = 0;
%     Bx = -30;
%     By = 10;
% 
    % Calculate the slope of both paths
    m1 = (By - Ay) / (Bx - Ax);
    m2 = -1 / m1;

    % Calculate the offset of both paths
    offset1 = Ay - (m1 * Ax);
    offset2 = Py - (m2 * Px);

    % Calculate intersection point
    xi = (offset1 - offset2) / (m2 - m1);
    yi = m2*xi + offset2;

    % Calculate image point
    x = (2 * xi) -Px;
    y = (2 * yi) -Py;
    imagePosition = [x y];

%     plot([Ax Bx], [Ay By]);
%     hold on
%     grid on
%     plot(xi, yi, 'x')
%     plot(Px, Py, 'o')
%     plot(x, y, 'd')
%     pause
%     xlim([-15 15])
%     ylim([-15 15])
%     return


end
