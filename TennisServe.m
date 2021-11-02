%% *TennisServe*
%% *definition*
classdef TennisServe
    %% *purpose*
    % to model the path of a tennis ball that has been served
    %% *description*
    % The user defines the initial position, velocity, and spin of the
    % tennis ball.  The units are pounds, feet, seconds
    % Tennis court coordinates
    %   x - ground, perpendicular to net
    %   y - ground, parallel to net
    %   z - height above ground
    % The path taken depends on 3 forces:
    %   1) gravity (constant 32 ft/s^2 in -z direction
    %   2) drag (proportional to v^2 and in -v-hat direction
    %   3) lift (proportional to cross product of spin and velocity)
    % The proportionality constant of the drag force is proportiaonl to the
    % atmosphere denisty, the cross-sectional area of the ball, the
    % coefficient of drag, and the inverse of mass:
    %                k_D = C_D*rho_atmosphere*area/(2*mass)
    % Spin of the ball causes it to lift. The constant has similar form:
    %                k_L = C_L*rho_atmosphere*area/(2*mass)
    % Spin will increase the margin of error when serving.  A top-spin
    % (+y-component) creates a downward force causing the ball to drop and
    % allowing the server to hit the ball at a greater initial angle.
    % Backspin is created by setting the y-component to a negative value.
    %% *references*
    % based on a model described in 
    % Danby, J. M. A. (1997) Computer Modeling: From Sports to Spaceflight 
    % ... From Order To Chaos, Willmann-Bell, Inc., Richmond, VA., p. 159.
    %% *joke of the day*
    % Q: What's an 'Open Set'?
    % A: Something that Pete Sampras usually wins at Forrest Hills
    % Q: What is |Rogers-x| < epsilon?
    % A: Mr. Roger's Neighborhood
    %% *history*
    %  When       Who    What
    %  ---------- ------ --------------------------------------------------
    %  2019/02/17 mnoah  original code
    properties(Constant)
        
        % acceleration due to gravity at Earth's surface
        g = 32; % (ft/s^2) 
        
        % parameters for the Runge-Kutta-Fehlberg method of order 4/5
        neq = 6; % number of equations
        a = [0 2/9 1/3 3/4 1 5/6];
        ch = [47/450  0  12/25 32/225 1/30 6/25];
        ct = [-1/150 0 3/100 -16/75 -1/20 6/25];
        b = [ ...
            [0   0 0 0 0]; ...
            [2/9 0 0 0 0]; ...
            [1/12 1/4 0 0 0]; ...
            [69/128 -243/128 135/64 0 0]; ...
            [-17/12 27/4 -27/5 16/15 0]; ...
            [65/432 -5/16 13/16 4/27 5/144]];
        
        % tennis court size constraints
        x_min = -39; % (ft) minimum distance from player to net
        x_max = 39;  % (ft) maximum distance from player to net
        y_min = -18; % (ft) half-width of tennis court and net
        y_max = 18;  % (ft) half-width of tennis court and net
        z_min = 0;   % (ft) the ground level
        
    end
    
    properties
        
        % about the tennis ball
        ball_radius =  1.28/12.0; % (ft) *2.54*1e-3; % (m)
        ball_mass = 0.13; % (lb)
        % coefficient of drag
        % (may not be constant for some flights, but treated
        % as such for this model)
        C_D = 0.4;   % (non-dimensional and unitless) [0 to 1]
        % coefficient of lift
        C_L = 0.107; % (non-dimensional and unitless) [0 to 1]
        % coefficient of restitution is
        % 1 = bounce height is same as max height of flight
        % for impact speed of 22 ft/s, restitution ~ 0.7
        % for impact speed of 81 ft/s, restitution ~ 0.5
        restitution = 0.5; % (non-dimensional and unitless) [0 to 1]
        % the grip_factor modifies the slippage relative to the surface
        % when a spinning ball collides with a surface       
        grip_factor = 0; % (non-dimensional and unitless) [0 to 1]
        
        % environmental parameters
        air_density = 0.075; % (pounds per cubic foot)
        
        % initial conditions of the serve
        % position
        x0 = -39;
        y0 = 8;
        z0 = 8;
        % velocity
        v0 = 150;       % (ft/s) initial ball speed
        theta0 = -4.7;  % (deg) angle of serve with respect to horizontal
        phi0 = 0.0;     % (deg) angle of serve with respect to x-axis
        % spin, in this model, is treated as a constant
        spin_x = 0; % (rad/s)
        spin_y = 0; % (rad/s)
        spin_z = 0; % (rad/s)
        
        % dynamic parameters of the simulation
        dt = 0.01; % (s) delta time (step-size) for propagation
        tol = 0.0000001; % tolerance difference for Runge-Kutta
                         % can change for different dynamics
        
    end
    
    properties(Dependent)
        ball_area; % (ft^2) cross-sectional ball area
        k_D; % (ft^-1) force constant due to drag
        k_L; % (ft^-1) force constant due to lift
        vx0; % (ft/s)
        vy0; % (ft/s)
        vz0; % (ft/s)
    end
    
    methods
        
        % methods to set the dependent parameters
        function [val] = get.ball_area(this)
            val = pi*this.ball_radius^2;
        end
        
        function [val] = get.k_D(this)
            val = this.C_D*this.air_density*this.ball_area/(2*this.ball_mass);
        end
        
        function [val] = get.k_L(this)
            val = this.C_L*this.air_density*this.ball_area/(2*this.ball_mass);
        end
        
        function [val] = get.vx0(this)
            val = this.v0*cosd(this.theta0)*cosd(this.phi0);
        end
        
        function [val] = get.vy0(this)
            val = this.v0*cosd(this.theta0)*sind(this.phi0);
        end
        
        function [val] = get.vz0(this)
            val = this.v0*sind(this.theta0);
        end
        
        % methods to set the user-defined parameters, with range checking
        function [this] = set.x0(this,inval)
            if (abs(inval) > 39)
                disp('[x0,y0,z0] (ft) are the coordinates of the tennis ball at time of serve.');
                error('[ERROR] abs(x) must be less than 39 ft, the length of half of the court');
            end
            this.x0 = inval;
        end
        
        function [this] = set.y0(this,inval)
            if (abs(inval) > 18)
                disp('[x0,y0,z0] (ft) are the coordinates of the tennis ball at time of serve.');
                error('[ERROR] abs(y) must be less than 18 ft, the half-width of the court');
            end
            this.y0 = inval;
        end
        
        function [this] = set.z0(this,inval)
            if (2 <= inval && inval <= 15)
                this.z0 = inval;
            else
                disp('[x0,y0,z0] (ft) are the coordinates of the tennis ball at time of serve.');
                error('[ERROR] tennis serve height allowed values: 2 ft to 15 ft');
            end
        end
        
        function [this] = set.v0(this,inval)
            if (0 <= inval && inval <= 250)
                this.v0 = inval;
            else
                disp('v0 (ft/s) is the initial speed of the tennis ball.');
                error('[ERROR] allowed initial speeds between 0 and 250 ft/s');
            end
        end
        
        function [this] = set.theta0(this,inval)
            if (-10 <= inval && inval <= 45)
                this.theta0 = inval;
            else
                disp('theta0 is the angle to the horizontal and phi0 is the angle with the x-axis.');
                error('[ERROR] allowed range -10 to 45 degrees wrt horizontal');
            end
        end
        
        function [this] = set.phi0(this,inval)
            if (-45 <= inval && inval <= 45)
                this.phi0 = inval;
            else
                disp('theta0 is the angle to the horizontal and phi0 is the angle with the x-axis.');
                error('[ERROR] allowed range -45 to 45 degrees wrt x-axis');
            end
        end
        
        function [this] = set.C_D(this,inval)
            if (0 <= inval && inval <= 1)
                this.C_D = inval;
            else
                disp('C_D is the non-dimensional and unitless coefficient of drag');
                error('[ERROR] allowed range 0 to 1');
            end
        end
        
        function [this] = set.C_L(this,inval)
            if (~exist('inval','var'))
                % for small spheres, the coefficient of lift is
                % roughly the ball radius in feet
                % (non-dimensional and unitless) [0 to 1] coefficient of lift
                estimated_value = this.ball_radius;
                inval = estimated_value;
            end
            if (0 <= inval && inval <= 1)
                this.C_L = inval;
            else
                disp('C_L is the non-dimensional and unitless coefficient of lift');
                error('[ERROR] allowed range 0 to 1');
            end
        end
        
        function [inputParm] = getInputs(this)
            
            inputParm.ball_radius = 1.28/12.0; 
            inputParm.ball_radius_unit = 'ft';
            inputParm.ball_mass = 0.13;  
            inputParm.ball_mass_unit = 'lb';
            inputParm.C_D = 0.4;
            inputParm.C_D_unit = 'unitless';
            inputParm.C_L = 0.107;
            inputParm.C_L_unit = 'unitless';
            inputParm.k_D = this.k_D;
            inputParm.k_D_unit = 'ft^-1';
            inputParm.k_L = this.k_L;
            inputParm.k_L_unit = 'ft^-1';
            inputParm.restitution = 0.5;
            inputParm.restitution_unit = 'unitless';
            inputParm.grip_factor = 0;
            inputParm.grip_factor_unit = 'unitless';
            inputParm.air_density = 0.075;
            inputParm.air_density_unit = 'lbs/ft^3';
            inputParm.initial_position = [this.x0, this.y0, this.z0];
            inputParm.initial_position_unit = 'ft';
            inputParm.initial_speed = this.v0;
            inputParm.initial_speed_unit = 'ft/s';
            inputParm.initial_velocity = [this.vx0, this.vy0, this.vz0];
            inputParm.initial_velocity_unit = 'ft/s';
            inputParm.initial_theta = this.theta0;
            inputParm.initial_theta_unit = 'deg';
            inputParm.initial_phi = this.phi0;
            inputParm.initial_phi_unit = 'deg';
            inputParm.initial_spin = [this.spin_x, this.spin_y, this.spin_z];
            inputParm.initial_spin_unit = 'rad/s';           

        end
        
        function [y] = Init(this)
            y = [this.x0 this.vx0 this.y0 this.vy0 this.z0 this.vz0];
        end
        
        function [z] = getAcceleration(this, y)
            Speed = sqrt(y(2)^2 + y(4)^2 + y(6)^2); % (ft/s)
            z(1) = y(2);
            z(2) = -this.k_D*Speed*y(2) + this.k_L*(this.spin_y*y(6) - this.spin_z*y(4));
            z(3) = y(4);
            z(4) = -this.k_D*Speed*y(4) + this.k_L*(this.spin_z*y(2) - this.spin_x*y(6));
            z(5) = y(6);
            z(6) = -this.k_D*Speed*y(6) + this.k_L*(this.spin_x*y(4) - this.spin_y*y(2)) - this.g;
        end
        
        function [Time_out,StepSize_out,y_out] = updateStateVector(this, Time, StepSize, y_in)
            
            Ttemp = Time;
            h = StepSize;
            f = zeros(6,6);
            f(1,:) = getAcceleration(this, y_in);
            ytemp = y_in;
            temax = 1e300;
            te = zeros(1,this.neq);
            
            % flagRun = true;
            while (temax >= this.tol)
                %fprintf(1,'%g %g\n',temax, this.tol);
                for k = 2:6
                    x = Ttemp + this.a(k)*h;
                    y = ytemp;
                    for n = 1:this.neq
                        for m = 1:k-1
                            y(n) = y(n) + h*this.b(k,m)*f(m,n);
                        end
                    end % n
                    f(k,:) = getAcceleration(this, y);
                end % k
                y = ytemp;
                % te = zeros(1,this.neq);
                for n = 1:this.neq
                    te(n) = 0;
                    for k = 1:6
                        y(n) = y(n) + h*this.ch(k)*f(k,n);
                        te(n) = te(n) + h*this.ct(k)*f(k,n);
                    end % k
                    te(n) = abs(te(n));
                end % n
                temax = this.tol/10.0;
                for n = 1:this.neq
                    if temax < te(n)
                        temax = te(n);
                    end
                end % n
                htemp = h;
                h = 0.9*h*exp(log(this.tol/temax)/5);
            end % repeat until temax < tol;
            y_out = y;
            Time_out = Ttemp + htemp;
            StepSize_out = h;
            
        end
        
        function [trajTable] = getTrajectory(this)
            % *outputs*
            % trajTable VariableNames:
            %  tB  - (s) time array
            %  xB  - (ft) x-position from net center, perpendicular to net
            %  yB  - (ft) y-position from net center, parallel to net
            %  zB  - (ft) altitude above the ground
            %  vxB - (ft/s) x-component speed
            %  vyB - (ft/s) y-component speed
            %  vzB - (ft/s) z-component speed
            % trajTable UserData parameters:
            %  flagBounce - true if the ball bounced
            %  flagNet - true if the ball hit the net
            %  flagOut - true if the ball went out of bounds
            y = Init(this);
            fprintf(1,'Start Pos=[%f,%f,%f] ft Vel=[%f,%f,%f] ft/s\n', ...
                y(1),y(3),y(5),y(2),y(4),y(6));
            Time = 0;
            StepNumber = 1;
            flagNet = false;
            flagOut = false;
            flagBounce = false;
            xB(1,1) = y(1); yB(1,1) = y(3); zB(1,1) = y(5); % ft
            tB(1,1) = Time;
            eps_bounce = 0.000000001;
            eps_trajectory = 0.0000001;
            eps_test = 0.00001;
            this.tol = eps_trajectory;
            
            xB(StepNumber,1)  = y(1);
            vxB(StepNumber,1) = y(2);
            yB(StepNumber,1)  = y(3);
            vyB(StepNumber,1) = y(4);
            zB(StepNumber,1)  = y(5);
            vzB(StepNumber,1) = y(6);
            tB(StepNumber,1) = Time; 
            
            while (~flagNet && ~flagOut && y(1) < this.x_max && StepNumber < 10000)
                
                StepNumber = StepNumber + 1;
                %fprintf(1,'%d\n',StepNumber);
                
                NextTime = Time + this.dt;
                while ((abs(NextTime-Time) > eps_test))
                    StepSize = NextTime - Time;
                    %fprintf(1,'%d %g\n',StepNumber,StepSize);
                    [Time,~,y] = updateStateVector(this, Time, StepSize, y);
                end

                xB(StepNumber,1)  = y(1);
                vxB(StepNumber,1) = y(2);
                yB(StepNumber,1)  = y(3);
                vyB(StepNumber,1) = y(4);
                zB(StepNumber,1)  = y(5);
                vzB(StepNumber,1) = y(6);
                current_time = tB(StepNumber-1,1) + StepSize;
                tB(StepNumber,1) = current_time;

                % kinematic solution
%  aaX = xB(StepNumber-1,1) + vxB(StepNumber-1,1)*StepSize;
%  aaY = yB(StepNumber-1,1) + vyB(StepNumber-1,1)*StepSize;
%  aaZ = zB(StepNumber-1,1) + vzB(StepNumber-1,1)*StepSize - 0.5*this.g*StepSize^2;
%  fprintf(1,'%d %f Pos=[%f,%f,%f] ft Vel=[%f,%f,%f] ft/s %f %f %f\n', ...
%    StepNumber,current_time,y(1),y(3),y(5),y(2),y(4),y(6),aaX,aaY,aaZ);
                            
                if ( (y(1)*xB(StepNumber-1,1) < 0) && (y(5) < 3) )
                    % ball is on server's side, but its height is below net
                    flagNet = true;
                    xB(StepNumber,1) = 0;
                    zB(StepNumber,1) = y(5);
                    StepNumber = StepNumber + 1;
                    Time = Time + 0.01; % mnoah
                    zB(StepNumber,1) = 0;
                    xB(StepNumber,1) = -y(2)*Time/10.0;
                    yB(StepNumber,1) = yB(StepNumber-1,1) + y(4)*Time/10;
                    
                    tB(StepNumber,1) = Time;
                    vxB(StepNumber,1) = 0;
                    vyB(StepNumber,1) = 0; 
                    vzB(StepNumber,1) = 0;
                elseif (abs(y(3)) > this.y_max)
                    % ball is outside the court (left-right)
                    flagOut = true;
                elseif (0 <= y(5) && y(5) < this.ball_radius)
                    % getting close to a bounce
                    this.tol = eps_bounce;
                elseif (y(5) < 0)
                    % the ball hit the ground, find the bounce velocity and
                    % position
                    this.tol = eps_bounce;
                    while (this.tol > abs(y(5)))
                        StepSize = -y(5)/y(6);
                        [Time,~,y] = updateStateVector(this, Time, StepSize, y);
                    end
                    y(2) = y(2) + this.grip_factor*this.spin_y;
                    y(4) = y(4) - this.grip_factor*this.spin_x;
                    y(6) = -this.restitution*y(6);

                    xB(StepNumber,1)  = y(1);
                    vxB(StepNumber,1) = y(2);
                    yB(StepNumber,1)  = y(3);
                    vyB(StepNumber,1) = y(4);
                    zB(StepNumber,1)  = 0.0;
                    vzB(StepNumber,1) = y(6);
                    tB(StepNumber,1) = tB(StepNumber-1,1) + StepSize;
                    this.tol = eps_trajectory;
                    flagBounce = true;
                end
                
            end % while ~flagNet ~flagOut and (xB < x_max)
            
            trajTable = table(tB,xB,yB,zB,vxB,vyB,vzB,'VariableNames', ...
                {'time','x','y','z','vx','vy','vz'});
            trajTable.Properties.VariableUnits = ...
                {'s','ft','ft','ft','ft/s','ft/s','ft/s'};
            trajTable.Properties.UserData.flagBounce = flagBounce;
            trajTable.Properties.UserData.flagOut = flagOut;
            trajTable.Properties.UserData.flagNet = flagNet;
            trajTable.Properties.UserData.inputParm = getInputs(this);
        end
        
    end
    
    
end

