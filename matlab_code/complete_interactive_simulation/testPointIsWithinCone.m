function result = testPointIsWithinCone(sourcePosition, boundary, receiverPosition)

    Px = sourcePosition(1);
    Py = sourcePosition(2);
    Rx = receiverPosition(1);
    Ry = receiverPosition(2);
    Ax = boundary(1);
    Ay = boundary(2);
    Bx = boundary(3);
    By = boundary(4);

    % Calculate cone's angle
    angle1 = atan( (Py - Ay) / (Px - Ax) );
    angle2 = atan( (Py - By) / (Px - Bx) );
    coneAngle = wrapToPi(angle2 - angle1);

    % Calculate cone's central vector (normalised)
    centralAngle = (angle1 + angle2) / 2;
    centralVector = [cos(centralAngle) sin(centralAngle)];

    % Calculate vector from cone's center to receiver's position
    sourceAngle = atan( (Ry - Py) / (Rx - Px) );
    queryVector = [cos(sourceAngle) sin(sourceAngle)];

    % Calculate the dot-product between vectors and take the arc-cos
    vectorsAngle = acos(dot(queryVector, centralVector));

    % Compare the angle to half the cone's angle to decide if it's within
    if vectorsAngle > coneAngle/2
        % It's outside
        result = 0;
    else
        % It's inside
        result = 1;
    end

%     plot(Px, Py, 'x', 'MarkerSize', 10, 'LineWidth', 3)
%     hold on
%     grid on
%     plot([boundary(1) boundary(3)], [boundary(2) boundary(4)], 'LineWidth', 2)
%     plot(Rx, Ry, 'o', 'MarkerSize', 10, 'LineWidth', 3)
%     plot([Px (5*centralVector(1))+Px], [Py (5*centralVector(2))+Py], '--')
%     plot([Px (5 * queryVector(1))+Px], [Py (5 * queryVector(2))+Py], '--')
%     plot([Px Ax], [Py Ay], ':')
%     plot([Px Bx], [Py By], ':')
%     if result
%         title('Is it within the cone: yes.');
%     else
%         title('Is it within the cone: no');
%     end
%     legend({'Source position', 'Boundary', 'Receiver position', 'Central vector', 'Query vector'}, 'Location', 'best')
    
end
