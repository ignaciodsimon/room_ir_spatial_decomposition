function boundaries = circular_room(radius, amount_segments)
    % Function to generate boudary matrices following a circumference built
    % with linear segments.
    %
    % Joe.

    circular_room_boundaries = zeros(amount_segments, 4);
    arc_angle = 360 / 180 * pi / amount_segments;
    for i = 1 : amount_segments
        circular_room_boundaries(i,:) = [[cos(arc_angle * i) sin(arc_angle * i)*radius] [cos(arc_angle * (i+1)) sin(arc_angle * (i+1))*radius]];
    end

    for boundary = circular_room_boundaries'
        plot([boundary(1) boundary(3)], [boundary(2) boundary(4)], 'Color', [.2 .8 .3])
        hold on
    end
end