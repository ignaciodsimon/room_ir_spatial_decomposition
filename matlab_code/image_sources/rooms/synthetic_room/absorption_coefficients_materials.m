function absorption_coefficients_materials()
    % Data from the book "Environmental and architectural acoustics"
    % by Z. Maekawa and P. Lord (from table A.2)
    %
    % NOTE: remember the original table lists absorption, not reflection
    % coefficients.
    %
    % Joe.

    % Soft urethane foam:
    porous = [ 125, 1-0.07
               250, 1-0.20
               500, 1-0.40
              1000, 1-0.55
              2000, 1-0.70
              4000, 1-0.70];

    % Brick, bare concrete surface
    wall = [ 125, 1-0.01
             250, 1-0.02
             500, 1-0.02
            1000, 1-0.02
            2000, 1-0.03
            4000, 1-0.04];

    % Wood-floor (parquet, or flooring on stud)
    wood = [ 125, 1-0.16
             250, 1-0.14
             500, 1-0.11
            1000, 1-0.08
            2000, 1-0.08
            4000, 1-0.07];

    % Window-glass (in wooden-frame)
    window = [ 125, 1-0.35
               250, 1-0.25
               500, 1-0.18
              1000, 1-0.12
              2000, 1-0.07
              4000, 1-0.04];

    % "2 mm steel flash door with air layer 45 mm", transmission loss (dB)
    % Converted to absorption coefficient following expression 5.26 on
    % page 162 of "Environmental and Architectural Acoustics", 1994.
    door = [ 125, 1-(10^(-25/10))
             250, 1-(10^(-30/10))
             500, 1-(10^(-34/10))
            1000, 1-(10^(-37/10))
            2000, 1-(10^(-36/10))
            4000, 1-(10^(-35/10))];

%     disp('Porous:');
%     disp(porous(:, 2)');
%     disp('Wall:')
%     disp(wall(:, 2)')
%     disp('Wood:')
%     disp(wood(:, 2)')
%     disp('Window:')
%     disp(window(:, 2)')
%     disp('Door:')
%     disp(door(:, 2)')
%     return

    % Create plot with all data
    figureHandler = figure();
    semilogx(porous(:,1), 1-porous(:,2), 'x-', 'LineWidth', 2)
    hold on
    semilogx(wall(:,1), 1-wall(:,2), 'x-', 'LineWidth', 2)
    semilogx(wood(:,1), 1-wood(:,2), 'x-', 'LineWidth', 2)
    semilogx(window(:,1), 1-window(:,2), 'x-', 'LineWidth', 2)
    semilogx(door(:,1), 1-door(:,2), 'x-', 'LineWidth', 2)
    legend({'Porous material', 'Concrete wall', 'Wood', 'Window', 'Door'}, 'Location', 'best')
    grid
    xlabel('Frequency [Hz]');
    ylabel('Absorption coefficient [.]');

    % Save plot to PDF
    set(gcf, 'PaperPosition', [-0.42 +0.05 9.3 3.6]);
    set(gcf, 'PaperSize', [8.1 3.5]);
    saveas(gcf, 'absorption_coefficients_materials.pdf', 'pdf');
    close(figureHandler);

end
