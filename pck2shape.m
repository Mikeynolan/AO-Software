function [angle0, angle1, angle2] = pck2shape(RA0, DEC0, W0)
    %
    % This will give values for an epoch of J2000. It could do a more
    % useful epoch, or that could happen elsewhere. Note that J2000 isn't a
    % legal date in shape, which requires integer seconds. Since it's J2000,
    % don't care about W1 and W2 
    
    %I am as literally as possible undoing shape2pck, which comes from
    %Steve Chesley's compute_W0
    
    
    % Equatorial vs. Ecliptic. Shape uses this value (IAU2006)
    eps = 84381.406/3600*pi/180;
    rot2eq = [1 0 0 ; 0 cos(eps) sin(eps); 0 -sin(eps) cos(eps)]';
    
    z_eq = [cosd(RA0)*cosd(DEC0), sind(RA0)*cosd(DEC0), sind(DEC0)];
    
    % Bennu equinox frame: X = equinox, Z = spin pole, Y = Z x X
    x_eqx_eq = cross([0 0 1], z_eq); x_eqx_eq = x_eqx_eq/norm(x_eqx_eq);
    y_eqx_eq = cross(z_eq, x_eqx_eq);
    z_eqx_eq = z_eq;
    rot_eq2eqx = [x_eqx_eq;y_eqx_eq;z_eqx_eq];
    
    %Bennu x-axis: Xz is 0 in both coordinate systems because it is by
    %definition the cross-product of the Zs
    
    x_bennu_eqx = [cosd(W0) sind(W0) 0];
    
    x_eq = rot_eq2eqx' * x_bennu_eqx';
    
    y_eq = cross(z_eq, x_eq);
    
    % Was rot2eq' in the other direction
    z_ec = z_eq * rot2eq;
    x_ec = x_eq' * rot2eq;
    y_ec = y_eq * rot2eq;
    
    r = [x_ec', y_ec', z_ec'];
    
    theta = acosd(r(3,3));
    psi = atan2d(r(3,1),r(3,2));
    phi = atan2d(r(1,3),-r(2,3));
    angle0 = phi;
    angle1=theta;
    angle2=psi;
end
