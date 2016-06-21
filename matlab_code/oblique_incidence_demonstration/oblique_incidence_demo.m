function oblique_incidence_demo()

    % Based on the formula on page Environmental and Architectural
    % acoustics, page 121, "sound absorption, materials and construction",
    % 1994.

    % Material impedance
    rn = 10;
    xn = 10;

    % Formula computation
    incidence_angles = [0 : 0.005 : pi/2 + 0.01];
    absorption = (4*rn*cos(incidence_angles)) ./ ((rn*cos(incidence_angles) + 1).^2 + (xn*cos(incidence_angles)).^2);

    % Results
    plot(incidence_angles / pi * 180, absorption, 'LineWidth', 1.5);
    grid
    xlim([0 91])
    ylim([0 1])
    xlabel('Incident angle rel. to normal incidence [degrees].')
    ylabel('Absorption coefficient [.]')
    title('Example of oblique incidence effect')

    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 0.15 9.4 3.25]);
    set(gcf, 'PaperSize', [8.5 3.5]);
    saveas(gcf, 'oblique_incidence_demo.pdf', 'pdf');
    close all
end
