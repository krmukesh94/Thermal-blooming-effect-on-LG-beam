%% Function of Initial Beam architecture of 44 Beam
% This will give the X and Y coordinate of centre of each laser source 
%%
function [p_x, p_y] =  Beam_Architecture1(N)
%% Default Values

N = 1;    % No of Laser 
r = 6e-2;   % Outer laser radius (machinical support)
R = 5e-2;     % Inner laser radius (Actual size oe laser beam)
l=6.5e-2;
origin_x = 0;
origin_y = 0;


%% Architecture:: The position of the beams
hold on;

    for i = 1:1;
        p_x(i) = origin_x;
        p_y(i) = origin_y;
        viscircles([p_x(i) p_y(i)] ,r);
        viscircles([p_x(i) p_y(i)] ,R, 'color' , 'blue');
        
    end

    
    
grid;
xlabel('Radial Distance (mm)'); ylabel('Radial Distance (mm)');
%title("Beam architecture")
xlim([-0.25 0.25])
ylim([-0.25 0.25])
