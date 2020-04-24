%% Read Forces.out file from AVL
% clear
% clc
%[FRED] Use this to post-process, you could make a similar file to read the
%output from the ST command

%INPUT: AVL forces.out file
%OUTPUT: forces .alpha .cl .cd .cm
% file = 'Trim_fold2';
function [Geo, Aero, Con] = readAVLForces(file)

fid = fopen(file);

while feof(fid)~=1
    line = fgetl(fid);
    
    
    % Inertial Data
    
    
    % %     Get the mass file and read that
    
    %     FlightData.Inertial.g = 9.81;           % Gravity Constant
    %     FlightData.Inertial.m = 2087;           % Aircraft Mass (kg)
    %     FlightData.Inertial.Ixx = 5066;         % Aircraft Moments of Inertia (kg.m^2)
    %     FlightData.Inertial.Iyy = 6578;         % Aircraft Moments of Inertia (kg.m^2)
    %     FlightData.Inertial.Izz = 10975;        % Aircraft Moments of Inertia (kg.m^2)
    %     FlightData.Inertial.Ixz = 203;          % Aircraft Moments of Inertia (kg.m^2)
    
    % Geometric Data
    if contains(line,'Sref')
        Geo.S = str2double(line(10:19));   % Platform Area (m^2)
    end
    if contains(line,'Cref')
        Geo.c = str2double(line(31:38));   % Chord Length (m)
    end
    if contains(line,'Bref')
        Geo.b = str2double(line(52:58));   % Wing Span (m)
    end
    
    % Aerodynamic Data (Reference CG:  % mac)
    %find Alpha
    if contains(line,'Alpha')
        Aero.Alpha = str2double(line(10:20));
    end
    %find CL
    if  contains(line,'CLtot')
        Aero.CL = str2double(line(10:19));
    end
    %find CD
    if  contains(line,'CDtot')
        Aero.CD = str2double(line(10:19));
    end
    %find Mach number
    if contains(line,'Mach')
        Aero.MN = str2double(line(10:19));
    end
    
    % Lift Coefficients
    if contains(line,'CLa')
        Aero.CLa = str2double(line(25:34));
    end
    if contains(line,'CLq')
        Aero.CLq = str2double(line(45:54));
    end
    if contains(line,'CLd3')
        Aero.CLde = str2double(line(65:74));
    end
    if  contains(line,'CLr')
        Aero.CLr = str2double(line(65:74));
    end
    
    % Side Force Coefficients
    if contains(line,'CYb')
        Aero.CYb = str2double(line(45:54));
    end
    if contains(line,'CYp')
        Aero.CYp = str2double(line(25:34));
    end
    if contains(line,'CYr')
        Aero.CYr = str2double(line(65:74));
    end
    if contains(line,'CYd2')
        Aero.CYda = str2double(line(45:54));
    end
    if contains(line,'CYd4')
        Aero.CYdr = str2double(line(85:94));
    end
    
    % M Moment Coefficients
    %find CM
    if  contains(line,'Cmtot')
        Aero.CM = str2double(line(32:41));
    end
    if  contains(line,'Cma')
        Aero.Cma = str2double(line(25:34));
    end
    if  contains(line,'Cmq')
        Aero.Cmq = str2double(line(45:54));
    end
    if  contains(line,'Cmd3')
        Aero.Cmde = str2double(line(65:74));
    end
    
    % N Moment Coefficients
    if  contains(line,'Cnb')
        Aero.Cnb = str2double(line(45:54));
    end
    if  contains(line,'Cnp')
        Aero.Cnp = str2double(line(25:34));
    end
    if  contains(line,'Cnr')
        Aero.Cnr = str2double(line(65:74));
    end
    if  contains(line,'Cnd2')
        Aero.Cnda = str2double(line(45:54));
    end
    if  contains(line,'Cnd4')
        Aero.Cndr = str2double(line(85:94));
    end
    
    % L Moment Coefficients
    if  contains(line,'Clb')
        Aero.Clb = str2double(line(45:54));
    end
    if  contains(line,'Clp')
        Aero.Clp = str2double(line(25:34));
    end
    if  contains(line,'Clr')
        Aero.Clr = str2double(line(65:74));
    end
    if  contains(line,'Cld2')
        Aero.Clda = str2double(line(45:54));
    end
    if  contains(line,'Cld4')
        Aero.Cldr = str2double(line(85:94));
    end
    
    % Control deflections 
    if  contains(line,'Flapd')
        Con.df = deg2rad(str2double(line(22:30)));
    end
    if  contains(line,'Ailerond')
        Con.da = deg2rad(str2double(line(22:30)));
    end
    if  contains(line,'Ruddervatored')
        Con.de = deg2rad(str2double(line(22:30)));
    end
    if  contains(line,'Ruddervatorrd')
        Con.dr = deg2rad(str2double(line(22:30)));
    end
    % Control surfaces Cd
    if  contains(line,'CDffd1')
        Con.Cddf = str2double(line(25:34));
    end
    if  contains(line,'CDffd2')
        Con.Cdda = str2double(line(45:54));
    end
    if  contains(line,'CDffd3')
        Con.Cdde = str2double(line(65:74));
    end
    if  contains(line,'CDffd4')
        Con.Cddr = str2double(line(85:94));
    end
end
end

