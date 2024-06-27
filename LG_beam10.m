%%
clc;
clear all; 
close all;

%% Beam Properties
N=1;   % No of Inputs
P = 1000;           % Power of each Beam(W) 
range = 0:500:5000;    % Slant range(n)
%range = 5000;
[px , py] = Beam_Architecture1(N);  % Call Function for Initial Beam architecture

%% Beam Parameters
lambda = 1064e-9;   % Wavelengtgh
W0 = 5e-2;            % Beam Waist radius(m)
rho = -0.2:0.005:0.2;  % Grid size
k = (2*pi)/lambda;  % Wave number
[X,Y] = meshgrid(rho, rho);	% Detector size

for n=1:length(range)
    %% Thermal Blooming
        
        alpha = 6.5*10^(-5);
        V = 2;             % Wind Velocity (m/s)
        n_T = -10^(-6);           %change in the refractive index with respect to temperature(1/K)
        n_0 = 1.000309;
        
        Cp = 10^3;      %specific heat of air at constant pressure of air(J/kg-K)
        rho_0 = 1.30175;   % ambient air density(kg/m^2)

    %% Laguerre Polynomiials, L_lm
    L_00 = @(x) 1;
    L_01 = @(x) 1 - x;
    L_02 = @(x) 1 - (2.*x) + ((x.^2)/2);
    L_03 = @(x) 8 .* (x.^3) - 12.*x;
    L_04 = @(x) 1 - (4.*x) + (3.*x.^3) - ((2/3).*x.^3) + (x.^4)./24;
    
    L_10 = @(x) 1;
    L_11 = @(x) 2 - x;
    L_12 = @(x) 3 - (3.*x) + ((x.^2)/2);
    L_13 = @(x) 4 - (9.*x) + (2.*x.^2) -(x.^4)./2;
    
    L_20 = @(x) 1;
    L_21 = @(x) 2 - x;
    L_22 = @(x) 3 - (3.*x) + ((x.^2)/2);
    
    
	I_combined = 0;
    for i = 1:N
        
        p_x(i) = px(i);
        p_y(i) = py(i);
      
        
        l=1;
        m=0;
        
        c1 = (1)/(2^(m+l)*factorial(l)*factorial(m)) * P/(pi*W0^2) ;
        I = c1 .* (exp((-2 .* ((X-p_x(i)).^2 + (Y-p_y(i)).^2))./W0.^2)) .* ((-2 .* ((X-p_x(i)).^2 + (Y-p_y(i)).^2))./W0.^2) .* (exp(-alpha.*range(n))) .* exp((((n_T.*alpha.*range(n).^2)./(2.*n_0.*rho_0.*Cp.*V))...
             .* (((-(1./(2 .*  W0.^5 .* ((X-p_x(i)).^2 + (Y-p_y(i)).^2)))) .* c1 .* exp(-((2 * ((X-p_x(i))^2 + (Y-p_y(i))^2))./W0.^2)) .* (4 .* W0 * (X-p_x(i)) * (W0.^2 * (5 * (X-p_x(i)).^2 + 7 * (Y-p_y(i)).^2) - 8 * ((X-p_x(i)).^4 + 3 * (X-p_x(i)).^2 * (Y-p_y(i))^2 + 2 * (Y-p_y(i))^4))...
             + (exp((2 .* (X-p_x(i))^2)./W0^2) * sqrt(2*pi) * (32 * (Y-p_y(i))^4 .* ((X-p_x(i)).^2 + (Y-p_y(i)).^2) + W0.^4 * ((X-p_x(i)).^2 + 3 * (Y-p_y(i))^2) - 4 .* W0.^2 .* (5 .* (X-p_x(i)).^2 .* (Y-p_y(i))^2 + 7 .* (Y-p_y(i))^4))) .* (1 + erf(sqrt(2).*(X-p_x(i))./W0)))))));
      
        I1 = c1 .* (exp((-2 .* ((X-p_x(i)).^2 + (Y-p_y(i)).^2))./W0.^2)) .* ((-2 .* ((X-p_x(i)).^2 + (Y-p_y(i)).^2))./W0.^2) .* (exp((-alpha.*range(n))) .* exp(((n_T.*alpha.*range(n).^2)./(2.*n_0.*rho_0.*Cp.*V))...
             .* (((-(1./(2 .*  W0.^5 .* ((X-p_x(i)).^2 + (Y-p_y(i)).^2)))) .* c1 .* exp(-((2 .* ((X-p_x(i)).^2 + (Y-p_y(i)).^2))./W0.^2)) .* (4 .* W0 .* (X-p_x(i)) .* (W0.^2 .* (5 .* (X-p_x(i)).^2 + 7 .* (Y-p_y(i)).^2) - 8 .* ((X-p_x(i)).^4 + 3 .* (X-p_x(i)).^2 .* (Y-p_y(i)).^2 + 2 .* (Y-p_y(i)).^4))...
             + (exp((2 .* (X-p_x(i)).^2)./W0.^2) .* sqrt(2.*pi) .* (32 .* (Y-p_y(i)).^4 .* ((X-p_x(i)).^2 + (Y-p_y(i)).^2) + W0.^4 .* ((X-p_x(i)).^2 + 3 .* (Y-p_y(i)).^2) - 4 .* W0.^2 .* (5 .* (X-p_x(i)).^2 .* (Y-p_y(i)).^2 + 7 .* (Y-p_y(i)).^4))) .* (1 + erf(sqrt(2).*(X-p_x(i))./W0)))))));
     
        I_combined= I_combined + I1;
        
    end
    I_combined = abs(I_combined)*1e-4;    % Convert in W / Cm2

    % Find Peak Intensity
    [maxI_com, rowp]  = max(I_combined);
    [Ipeak_com, Ipos] = max(maxI_com);
    peak_intensity = Ipeak_com
        
%% spot size
% spot size in the atmosphere

I_spot = @(x,y) ( c1 .* (exp((-2 .* (x.^2 + y.^2))./W0.^2)) .* ((-2 .* (x.^2 + y.^2))./W0.^2) .* (exp(-alpha.*range(n))) .* exp((((n_T.*alpha.*range(n).^2)./(2.*n_0.*rho_0.*Cp.*V))...
             .* (((-(1./(2 .*  W0.^5 .* (x.^2 + y.^2)))) .* c1 .* exp(-((2 .* (x.^2 + y.^2))./W0.^2)) .* (4 .* W0 .* x .* (W0.^2 .* (5 .* x.^2 + 7 .* y.^2) - 8 .* (x.^4 + 3 .* x.^2 .* y.^2 + 2 .* y.^4))...
             + (exp((2 .* x.^2)./W0.^2) .* sqrt(2.*pi) .* (32 .* y.^4 .* (x.^2 + y.^2) + W0.^4 .* (x.^2 + 3 .* y.^2) - 4 .* W0.^2 .* (5 .* x.^2 .* y.^2 + 7 .* y.^4))) .* (1 + erf(sqrt(2).*x./W0)))))))).^2;
     
I_spot_x = @(x,y) x.*(c1 .* (exp((-2 .* (x.^2 + y.^2))./W0.^2)) .* ((-2 .* (x.^2 + y.^2))./W0.^2) .* (exp(-alpha.*range(n))) .* exp((((n_T.*alpha.*range(n).^2)./(2.*n_0.*rho_0.*Cp.*V))...
             .* (((-(1./(2 .*  W0.^5 .* (x.^2 + y.^2)))) .* c1 .* exp(-((2 .* (x.^2 + y.^2))./W0.^2)) .* (4 .* W0 .* x .* (W0.^2 .* (5 .* x.^2 + 7 .* y.^2) - 8 .* (x.^4 + 3 .* x.^2 .* y.^2 + 2 .* y.^4))...
             + (exp((2 .* x.^2)./W0.^2) .* sqrt(2.*pi) .* (32 .* y.^4 .* (x.^2 + y.^2) + W0.^4 .* (x.^2 + 3 .* y.^2) - 4 .* W0.^2 .* (5 .* x.^2 .* y.^2 + 7 .* y.^4))) .* (1 + erf(sqrt(2).*x./W0)))))))).^2;
     
I_spot_y = @(x,y) y.*(c1 .* (exp((-2 .* (x.^2 + y.^2))./W0.^2)) .* ((-2 .* (x.^2 + y.^2))./W0.^2) .* (exp(-alpha.*range(n))) .* exp((((n_T.*alpha.*range(n).^2)./(2.*n_0.*rho_0.*Cp.*V))...
             .* (((-(1./(2 .*  W0.^5 .* (x.^2 + y.^2)))) .* c1 .* exp(-((2 .* (x.^2 + y.^2))./W0.^2)) .* (4 .* W0 .* x .* (W0.^2 .* (5 .* x.^2 + 7 .* y.^2) - 8 .* (x.^4 + 3 .* x.^2 .* y.^2 + 2 .* y.^4))...
             + (exp((2 .* x.^2)./W0.^2) .* sqrt(2.*pi) .* (32 .* y.^4 .* (x.^2 + y.^2) + W0.^4 .* (x.^2 + 3 .* y.^2) - 4 .* W0.^2 .* (5 .* x.^2 .* y.^2 + 7 .* y.^4))) .* (1 + erf(sqrt(2).*x./W0)))))))).^2;
     
%Shift in centriod
jx = integral2(I_spot_x,-inf, inf,-inf, inf)./integral2(I_spot,-inf, inf,-inf, inf)
jy = integral2(I_spot_y,-inf, inf,-inf, inf)./integral2(I_spot,-inf, inf,-inf, inf)

I_spot_Wx = @(x,y) (x-jx).^2 .*(c1 .* (exp((-2 .* (x.^2 + y.^2))./W0.^2)) .* ((-2 .* (x.^2 + y.^2))./W0.^2) .* (exp(-alpha.*range(n))) .* exp((((n_T.*alpha.*range(n).^2)./(2.*n_0.*rho_0.*Cp.*V))...
             .* (((-(1./(2 .*  W0.^5 .* (x.^2 + y.^2)))) .* c1 .* exp(-((2 .* (x.^2 + y.^2))./W0.^2)) .* (4 .* W0 .* x .* (W0.^2 .* (5 .* x.^2 + 7 .* y.^2) - 8 .* (x.^4 + 3 .* x.^2 .* y.^2 + 2 .* y.^4))...
             + (exp((2 .* x.^2)./W0.^2) .* sqrt(2.*pi) .* (32 .* y.^4 .* (x.^2 + y.^2) + W0.^4 .* (x.^2 + 3 .* y.^2) - 4 .* W0.^2 .* (5 .* x.^2 .* y.^2 + 7 .* y.^4))) .* (1 + erf(sqrt(2).*x./W0)))))))).^2;
     
I_spot_Wy = @(x,y) (y-jy).^2 .*(c1 .* (exp((-2 .* (x.^2 + y.^2))./W0.^2)) .* ((-2 .* (x.^2 + y.^2))./W0.^2) .* (exp(-alpha.*range(n))) .* exp((((n_T.*alpha.*range(n).^2)./(2.*n_0.*rho_0.*Cp.*V))...
             .* (((-(1./(2 .*  W0.^5 .* (x.^2 + y.^2)))) .* c1 .* exp(-((2 .* (x.^2 + y.^2))./W0.^2)) .* (4 .* W0 .* x .* (W0.^2 .* (5 .* x.^2 + 7 .* y.^2) - 8 .* (x.^4 + 3 .* x.^2 .* y.^2 + 2 .* y.^4))...
             + (exp((2 .* x.^2)./W0.^2) .* sqrt(2.*pi) .* (32 .* y.^4 .* (x.^2 + y.^2) + W0.^4 .* (x.^2 + 3 .* y.^2) - 4 .* W0.^2 .* (5 .* x.^2 .* y.^2 + 7 .* y.^4))) .* (1 + erf(sqrt(2).*x./W0)))))))).^2;
     
Wx = sqrt(4*integral2(I_spot_Wx,-inf, inf,-inf, inf)/integral2(I_spot,-inf, inf,-inf, inf))
Wy = sqrt(4*integral2(I_spot_Wy,-inf, inf,-inf, inf)/integral2(I_spot,-inf, inf,-inf, inf))

% spot size in free space 
I_free = @(x,y) (c1 .* (exp((-2 .* (x.^2 + y.^2))./W0.^2)) .* ((-2 .* (x.^2 + y.^2))./W0.^2)).^2;
I_free_x = @(x,y) x.*(c1 .* (exp((-2 .* (x.^2 + y.^2))./W0.^2)) .* ((-2 .* (x.^2 + y.^2))./W0.^2)).^2;
I_free_y = @(x,y) y.*(c1 .* (exp((-2 .* (x.^2 + y.^2))./W0.^2)) .* ((-2 .* (x.^2 + y.^2))./W0.^2)).^2;

jx_free = integral2(I_free_x,-inf, inf,-inf, inf)/integral2(I_free,-inf, inf,-inf, inf)
jy_free = integral2(I_free_y,-inf, inf,-inf, inf)/integral2(I_free,-inf, inf,-inf, inf)

I_free_Wx = @(x,y) (x-jx_free).^2 .*(c1 .* (exp((-2 .* (x.^2 + y.^2))./W0.^2)) .* ((-2 .* (x.^2 + y.^2))./W0.^2)).^2;
I_free_Wy = @(x,y) (y-jy_free).^2 .*(c1 .* (exp((-2 .* (x.^2 + y.^2))./W0.^2)) .* ((-2 .* (x.^2 + y.^2))./W0.^2)).^2;

Wx_free = sqrt(4*integral2(I_free_Wx,-inf, inf,-inf, inf)/integral2(I_free,-inf, inf,-inf, inf))
Wy_free = sqrt(4*integral2(I_free_Wy,-inf, inf,-inf, inf)/integral2(I_free,-inf, inf,-inf, inf))


%% Plot the output

    if(range(n) == 0 || range(n) == 100 || range(n) == 3000 || range(n) == 5000 || range(n) == 10000)
        figs(1)=figure;
        surf(X*100,Y*100,I_combined);
        xlabel('Radial Distance (cm)','fontweight','bold','fontsize',12); zlabel('Intensity (W/cm^2)','fontweight','bold','fontsize',12);
        grid; hcb=colorbar; view(3);
        colorTitleHandle = get(hcb,'Title');
            titleString = '{\it I(W/cm^2)}';
            set(colorTitleHandle ,'String',titleString);
        ax = gca; 
        set(gca,'FontWeight','bold');
        
        figs(2)=figure;
        image([-max(rho)*100 max(rho)*100], [-max(rho)*100 max(rho)*100], I_combined, 'CDataMapping','scaled')
        xlabel('X (cm)','fontweight','bold','fontsize',18); ylabel('Y (cm)','fontweight','bold','fontsize',18);
        hcb=colorbar;
            colorTitleHandle = get(hcb,'Title');
            titleString = '{\it I(W/cm^2)}';
            set(colorTitleHandle ,'String',titleString);
        ax = gca; 
        set(gca,'FontWeight','bold');
        
        figs(3)=figure;
        contour(X*100,Y*100,I_combined);
        xlabel('X (cm)','fontweight','bold','fontsize',18); ylabel('Y (cm)','fontweight','bold','fontsize',18);
        grid;hcb=colorbar;
        colorTitleHandle = get(hcb,'Title');
            titleString = '{\it I(W/cm^2)}';
            set(colorTitleHandle ,'String',titleString);
        ax = gca; 
        set(gca,'FontWeight','bold');
        %print figures to files
         for i_fig=1:length(figs)
            print(figs(i_fig), '-dpng', sprintf('LG10_[z=%dkm].png', N, range(n)/1000, i_fig), '-r300');
        end

   end
        %     dump values to excel file
        T = table(N, P, range(n), 10, Ipeak_com, jx, jy, Wx, Wy, Wx_free, Wy_free, Wx/Wx_free, Wy/Wy_free, Wy/Wx);
        filename = sprintf('Focused.xlsx');
        [d1,d2, existingData] = xlsread('Focused.xlsx');
        numberOfRow = size(existingData,1);     
        writetable(T,filename,'Sheet',1,'range',strcat('A',num2str(numberOfRow+1)),'WriteVariableNames',false); 
end



