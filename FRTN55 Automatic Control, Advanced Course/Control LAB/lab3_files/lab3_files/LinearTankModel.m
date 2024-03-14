classdef LinearTankModel
    %%LinearTankModel   Construct a linear model for the tank process
    %
    %   Construction:
    %       lt = LinearTankModel creates and empy LinearTankModel object
    %
    %   You set the properties
    %       TankArea = [A1 A2 A3 A4]
    %       OutletArea = [a1 a2 a3 a4]
    %       PumpConstants = [k1 k2]
    %       MeasurementConstant = [kc]
    %       Gammas = [gamma1, gamma2]
    %       Gravity, default is 981
    %       LinPointLowerTanks = [h1_0 h2_0]
    %       SampleTime
    %       PumpMax = [M1, M2]
    %       PumpMin = [m1, m2]
    %       TankHeights = [h1, h2, h3, h4]
    %       RateMin = [m m]
    %       RateMax = [M M]
    %
    %   You can then access the dependent properties which are
    %   automatically calculated
    %
    %   lt.T = [T1 T2 T3 T4]
    %   lt.CtSys is a continuous-time ss object corresponding to the
    %   specified properties
    %   lt.DtSys is a first-order hold discrete-time ss object
    %   corresponding to the specified properties.
    properties
        % Cross section of of each tank in cm^2
        % Set this property as a bj.LinPointLowerTanks(2)1x4 vector
        TankArea = 4.9*ones(1,4);
        % Cross section of each outlet in cm^2
        % Set this property as a 1x4 vector
        OutletArea = 0.03*ones(1,4);
        % Pump constants, in cm^3 V^{-1} s^{-1}
        % set as a 1x2 arraybj.LinPointLowerTanks(2)
        PumpConstants = 1.6*ones(1,2);
        % Measurement constant in V
        % Set this property as a scalar
        MeasurementConstant = 0.5;
        % The gamma values of the pumps
        % Set this value as a 1x2 vector
        Gammas = [.7, .7];
        % Gravity constant, defaults to 981 cm/s^2
        Gravity = 981
        % Linearization point for the lower tanks,
        % Set this property as a 1x2 vector
        LinPointLowerTanks = [10, 10];
        % SampleTime for discretization
        % Set this property as a scalar
        SampleTime = 1;
        % Maximum values for the pumps
        % set as a 1x2 array
        estimData = 0;
        useEstimData = 0;
    end
    properties (Dependent)
        T
        CtSys
        DtSys
        % Operating point on the form of [h1,h2,h3,h4,u1,u2]
        LinPoint
    end
    methods
        function LinPoint = get.LinPoint(obj)
            if obj.useEstimData == 0;
                % Rename variables for brevity
                h1 = obj.LinPointLowerTanks(1); h2 = obj.LinPointLowerTanks(2);
                g1 = obj.Gammas(1); g2 = obj.Gammas(2);
                k1 = obj.PumpConstants(1); k2 = obj.PumpConstants(2);
                a = obj.OutletArea;
                g = obj.Gravity;
                
                % Calculate the stationary points of actuator signal u and
                % lower tanks h3, h4
                
                u1 = sqrt(2*g)*((g2-1)*a(2)*sqrt(h2)+g2*a(1)*sqrt(h1))/(k1*(g1+g2-1)); % [V]
                u2 = sqrt(2*g)*((g1-1)*a(1)*sqrt(h1)+g1*a(2)*sqrt(h2))/(k2*(g1+g2-1)); % [V]
                h3 = (k2*u2*(1-g2)/a(3))^2/(2*g); % [cm]
                h4 = (k1*u1*(1-g1)/a(4))^2/(2*g); % [cm]
                LinPoint = [h1 h2 h3 h4 u1 u2];
            else
                h1 = obj.LinPointLowerTanks(1); h2 = obj.LinPointLowerTanks(2);
                B = obj.estimData.B;
                B1 = [B(1,1) B(3,2); B(4,1) B(2,2)];
                abyA = obj.estimData.abyA;
                u0 = B1\...
                    [abyA(1)*sqrt(2*obj.Gravity*h1)
                    abyA(2)*sqrt(2*obj.Gravity*h2)];
                h3 = (B(3,2)*u0(2)/abyA(3))^2/2/obj.Gravity;
                h4 = (B(4,1)*u0(1)/abyA(4))^2/2/obj.Gravity;
                LinPoint = [h1 h2 h3 h4 u0(1) u0(2)];
            end
        end
        function T = get.T(obj)
            if obj.useEstimData == 0;
                T = obj.TankArea./obj.OutletArea.*...
                    sqrt(2*obj.LinPoint(1:4)/obj.Gravity);
            else
                T = 1./obj.estimData.abyA.*...
                    sqrt(2*obj.LinPoint(1:4)/obj.Gravity);
            end
        end
        function CtSys = get.CtSys(obj)
            T = obj.T;
            a = obj.TankArea;
            g = obj.Gammas;
            k = obj.PumpConstants;
            A=[-1/T(1) 0 a(3)/(a(1)*T(3)) 0;
                0 -1/T(2) 0 a(4)/(a(2)*T(4));
                0 0 -1/T(3) 0;
                0 0 0 -1/T(4)];
            if obj.useEstimData == 0
            B=[g(1)*k(1)/a(1) 0;
                0 g(2)*k(2)/a(2);
                0 (1-g(2))*k(2)/a(3);
                (1-g(1))*k(1)/a(4) 0];
            else % in case the parameters are estimated
                B = obj.estimData.B;
            end
            C = obj.MeasurementConstant*eye(4);
            
            D=zeros(4,2);
            CtSys = ss(A, B, C, D);
            
            % Specify inputs, outputs and states
            CtSys.InputName = {'u_1', 'u_2'};
            CtSys.OutputName = {'y_1', 'y_2', 'y_3', 'y_4'};
            CtSys.StateName = {'h_1', 'h_2', 'h_3', 'h_4'};
            
        end
        
        function DtSys = get.DtSys(obj)
            DtSys = c2d(obj.CtSys, obj.SampleTime);
        end
    end
end
