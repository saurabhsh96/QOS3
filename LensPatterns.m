%Electric far field of in the infinite media
%clear;
function [EthInL, EphInL, PradFeed, rho, phi_s, JsX, JsY, JsZ, ...
    Axx, Ayy, EFLx, EFLy, EFLz, Prad, etaRad] = LensPatterns(freq, er, n, Dl, Th0, jf, ...
    dtho, dpho, theta_obs, phi_obs)
    %% Input definations

    %Frequency
    %freq = 100e9;

    %Material characteristics
    eps_0 = 8.854187817e-12;
    mu_0 = 1.2566370614e-6;
    %er = 11.9;

    %Speed of light in vaccum
    c = 3e8;

    %Propagation velocity in the medium
    vp = c/sqrt(er);

    %Impedance (in Ohm)
    zeta0 = round(sqrt(mu_0/(eps_0)));
    zetad = round(sqrt(mu_0/(eps_0*er)));

    %Wavelenth
    lam = c/freq;

    %Prop constnat
    k0 = 2*pi/lam;

    %Prop constnat in the medium
    kd = k0*sqrt(er);

    %Order (What exactly is this?)
    %n = 3;

    %% Field calculations
    %Defining the meshgrid (for Q1 it is entire half-space)
    drad = pi/180;
    [th, phi] = meshgrid(eps:drad:pi/2, eps:drad:2*pi);

    %Calling the Electric field function: Field in Media
    r_obs = 1000*lam;
    [Eth, Ephi] = FieldMedia(n, kd, th, phi, r_obs);
    
    %%%%%%%%%%%%
    EthInL = Eth;
    EphInL = Ephi;
    
    Emag = sqrt(abs(Eth).^2 + abs(Ephi).^2);
    Emax = max(max(Emag));

    %Calculate field intensity
    IntensityFeed = r_obs^2.*abs(Emag).^2./(2*zetad);

    %Calculate prad of feed 
    PradFeed = (sum(IntensityFeed.*sin(th), 'all'))*drad*drad; %dth = drad; dph = drad;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %InLInten = IntensityFeed;
    %InLPrad = PradFeed;
    
    %Plotting electric far field in the medium
%    figure(1);
% %    Phi = 0
%     plot([-th(1,91:-1:1)./drad, th(1,1:91)./drad], [mag2db(Emag(1,91:-1:1)/Emax), mag2db(Emag(1,1:91)/Emax)], 'LineWidth', 1.5);
%     title('Normalized Electric Field in \phi = 0 plane');
%     xlabel('\theta [Degrees]');
%     ylim([-60, 0]);
%     ylabel('Normalized Electric field [dBV]');
% 
%     figure(2);
%     % Phi = 0
%     plot([-th(1,91:-1:1)./drad, th(1,1:91)./drad], [mag2db(Emag(91,91:-1:1)/Emax), mag2db(Emag(91,1:91)/Emax)], 'LineWidth', 1.5);
%     title('Normalized Electric Field in \phi = 90 plane');
%     xlabel('\theta [Degrees]');
%     ylim([-60, 0]);
%     ylabel('Normalized Electric field [dBV]');

    %% Q2: Aperture Current distribution

    %Diameter of the lens
%    Dl = 6*lam;

    %Critical angle ASK
    %Th0 = 73.1489*drad;
    %Th0 = 56.2978*drad;
    %Th0 = 0.6658;
%    Th0 = 55*drad;

    %Eccentricity of ellipse
    e = 1/sqrt(er);
    %Th0 = asin(e);

    %Distance Variables
    r_min = Dl/2/sin(Th0);
    al = r_min*(1-e*cos(Th0))/(1-e^2);

    %Ellipse Parameters
    cent = al*e;
    b = sqrt(al^2 - cent^2);

    %Surface Parameterization
    dRho = 0.00001;
    dPhi_s = drad;
    rho_var = eps:dRho:Dl/2; %Send
    [rho, phi_s] = meshgrid(rho_var, eps:drad:2*pi); %Send
    z = sqrt((al.^2).*(1-(rho_var./b).^2))+cent;
    th_req_arr = atan(rho_var./z);

    %Transmission coefficient
    %Required meshgrid
    [th_req, phi_req] = meshgrid(th_req_arr, eps:drad:2*pi); 
    cos_in = (1 - e.*cos(th_req))./(sqrt(1+e^2-(2*e).*cos(th_req)));
    cos_tx = (cos(th_req)-e)./(sqrt(1+e^2-(2*e).*cos(th_req)));

    %Gamma calculations
    % gamma_par = (zeta0.*cos_in - zetad.*cos_tx)./(zeta0.*cos_in + zetad.*cos_tx);
    % figure();
    % plot(th_req(1,:)./drad, gamma_par(1,:));

    %Perpendicular Tx Coeff
    tau_perp = ((2*zeta0).*cos_in)./(zeta0.*cos_in + zetad.*cos_tx);
    tau_par = ((2*zeta0).*cos_in)./(zeta0.*cos_tx + zetad.*cos_in);
    Sr = sqrt((cos_tx./cos_in).*(e*cos(th_req)-1)./(e-cos(th_req)));

    %Field considering different phases
    [EthL, EphiL] = FieldMedia(n, kd, th_req, phi_req, r_obs);
    EmagL = sqrt(abs(EthL).^2 + abs(EphiL).^2);
    EmaxL = max(max(EmagL));

    %Remove initial phase
    mult_term = exp(-1j*kd*r_obs)./r_obs;
    %Add constant phase
    r = al.*(1-e^2)./(1-e.*cos(th_req));
    h = al.*(1+e); 
%    r_dash = h - r.*cos(th_req);
    %Should it be h in the constant term or r technically, the ray travels r+r'?
    %ASK this question, this effects the directivity and gain
    const_term = (Sr./(r)).*exp(-1j*kd*al*(1+e));
    Eth_wp = (EthL./mult_term).*const_term;
    Ephi_wp = (EphiL./mult_term).*const_term;
    EmagL_wp = sqrt(abs(Eth_wp).^2 + abs(Ephi_wp).^2);
    EmaxL_wp = max(max(EmagL_wp));

    %Current Distribution %Send JsX, JxY, JxZ
    JsX = (-2/zeta0).*(tau_par.*Eth_wp.*cos(phi_req) - tau_perp.*Ephi_wp.*sin(phi_req));
    JsY = (-2/zeta0).*(tau_par.*Eth_wp.*sin(phi_req) + tau_perp.*Ephi_wp.*cos(phi_req));
    JsZ = zeros(size(rho));

    %Conversion to X and Y send Axx, Ayy
    Axx = (rho).*cos(phi_s);
    Ayy = rho.*sin(phi_s);

    %Plotting Currents
%     figure(3);
%     surface(Axx, Ayy, mag2db(abs(JsX)./max(max(abs(JsX)))), 'linestyle','none');
%     title('Normalized Current Distribution at Equivalent Aperture, Jx');
%     xlabel('x [m]');
%     ylabel('y [m]');
% 
%     figure(4);
%     surface(Axx, Ayy, mag2db(abs(JsY)./max(max(abs(JsY)))), 'linestyle','none');
%     title('Normalized Current Distribution at Equivalent Aperture, Jy');
%     xlabel('x [m]');
%     ylabel('y [m]');

    %% Far field of the distribution

    %Meshgrid for observation
    %dtho = drad/8;
    %dpho = pi/180;
    %[theta_obs, phi_obs] = meshgrid(eps:dtho:pi/2-dtho, eps:dpho:2*pi); %send

    %JFT of the current, using function from last assignment
    [JFTx, JFTy, JFTz] = ApertureJFT(JsX, JsY, JsZ, dRho, dPhi_s, k0, rho, ...
        phi_s, theta_obs, phi_obs);

    JFTa = zeros([3, size(JFTx)]);
    JFTa(1,:,:)=JFTx;
    JFTa(2,:,:)=JFTy;
    JFTa(3,:,:)=JFTz;

    [EFLx, EFLy, EFLz] = FFRef(freq, 1, r_obs, theta_obs, phi_obs, JFTa);
    %Twice the current is considered and sent is again twice hence normalize it
    %by 2
    EFLx = EFLx./2; %send
    EFLy = EFLy./2; %send
    EFLz = EFLz./2; %send
    EFLMag = sqrt(abs(EFLx).^2 + abs(EFLy).^2 + abs(EFLz).^2); 
    EFLMagMax = max(max(EFLMag));

    %Field by unform aperture
%     jf = [1, 0, 0]';
    [EFxu, EFyu, EFzu] = FFRefFeed(freq, 1, jf, Dl/2, r_obs, theta_obs, phi_obs);
    EFRMagu = sqrt(abs(EFxu).^2 + abs(EFyu).^2 + abs(EFzu).^2); 
    EFRMagMaxu = max(max(EFRMagu));

%     figure(5);
%     plot([-theta_obs(1,720:-1:1)*180/pi, theta_obs(1,:)*180/pi], [mag2db(EFLMag(91,720:-1:1)/EFLMagMax), mag2db(EFLMag(91,:)/EFLMagMax)]); 
%     %hold on;
%     %plot([-theta_obs(1,720:-1:1)*180/pi, theta_obs(1,:)*180/pi], [mag2db(EFRMagu(91,720:-1:1)/EFRMagMaxu), mag2db(EFRMagu(91,:)/EFRMagMaxu)], '--');
%     title('Normalized E-field at plane \phi = 0');
%     ylabel('Normalized E-field [dBV]');
%     xlabel('Observation angle, \theta (in deg)');
%     ylim([-45, 0]);
%     %legend('Lens Aperture', 'Uniform Aperture');
% 
%     figure(6);
%     plot([-theta_obs(1,720:-1:1)*180/pi, theta_obs(1,:)*180/pi], [mag2db(EFLMag(91,720:-1:1)/EFLMagMax), mag2db(EFLMag(91,:)/EFLMagMax)]); 
%     %hold on;
%     %plot([-theta_obs(1,720:-1:1)*180/pi, theta_obs(1,:)*180/pi], [mag2db(EFRMagu(91,720:-1:1)/EFRMagMaxu), mag2db(EFRMagu(91,:)/EFRMagMaxu)], '--');
%     title('Normalized E-field at plane \phi = 90');
%     ylabel('Normalized E-field [dBV]');
%     xlabel('Observation angle, \theta (in deg)');
%     ylim([-35, 0]);
%     %legend('Lens Aperture', 'Uniform Aperture');

    %% Calculation of Directivity

    %Calculate field intensity
    Intensity = r_obs^2.*abs(EFLMag).^2./(2*zeta0);

    %Calculate prad
    Prad = (sum(Intensity.*sin(theta_obs), 'all'))*dtho*dpho; %dth = drad; dph = pi/4

    %Calculate directivity
    Dir = (4*pi).*Intensity/Prad;

    %Efficiency of lens
    etaRad = Prad/PradFeed;

    %Gain of the lens
    gain = Dir.*etaRad;

%     figure(7);
%     plot(theta_obs(1,:)./drad, pow2db(Dir(1,:)), '--'); hold on;
%     plot(theta_obs(1,:)./drad, pow2db(gain(1,:)));
%     title('Directivity and Gain at \phi = 0');
%     xlabel('\theta (Deg)');
%     ylabel('Gain [dBi] and Directivity [dBi]');
%     legend('D(\theta, 0)', 'G(\theta, 0)');
%     ylim([-30, 40]);

    %% Phase center evaluation
    dz = 0.0001;
    dth = drad;
    delZ = -3*lam:0.5*lam:3*lam;
    %phase = exp(1j*k0*cos(theta_obs).*delZ);

    %Parametric analysis over delZ
    %figNum = 8;
%    figure;
    for value = delZ
        EFLxn = EFLx.*exp(1j*k0*cos(theta_obs).*value);
        EFLyn = EFLy.*exp(1j*k0*cos(theta_obs).*value);
        phaseX = angle(EFLxn);
        phaseY = angle(EFLyn);
        %figure(figNum);
%         name = ['\Delta Z = ', num2str(value)];
%         plot(theta_obs(1,1:50)./drad, unwrap(angle(EFLxn(1,1:50))-angle(EFLxn(1,1)))./drad, 'DisplayName', name, 'LineWidth', 1.5); 
%         hold on;
%         %plot(theta_obs(1,1:50)./drad, unwrap(angle(EFLyn(1,1:50))-angle(EFLyn(1,1)))./drad);
%         %hold on;    
%         title('Phase of Ex vs. \Theta with varying \Delta Z');
%         xlabel('\theta [deg]');
%         ylabel('Phase angle [deg]');
        %legend('\Delta Z' + value);
        %figNum = figNum+1;
        %hold off;
    end
%    hold off;
%    legend show
end