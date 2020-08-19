% Code created by Loïc Marrec

function pfix = Gillespie_fct(Nit, n, theta, gA, XA_i, fB, gB, XB_i, K, t0)

    XB_list = NaN(1, Nit); 

    for i = 1 : Nit

        % Initialization

        w = 0;
        XA = XA_i;              % Initialization of the number of A individuals
        XB = XB_i;              % Initialization of the number of B individuals
        t = 0;                 % Initialization of time
        cumul = zeros(4, 1);    % To build the sampling tower

        while (XA ~= 0 && XB ~= 0) || w == 0
            
            % Compute the fitness of W microbes

            fA = sigm(t, theta, n);
            
            % Compute the transition rates

            T_A_rep = fA*(1-(XA+XB)/K)*XA;   
            T_A_death = gA*XA;
            T_B_rep = fB*(1-(XA+XB)/K)*XB;
            T_B_death = gB*XB;
            T = T_A_rep+T_A_death+T_B_rep+T_B_death;

            % Compute tau and then the new time

            r1 = rand;
            tau = 1/T*log(1/r1);
            t = t+tau;
            
            if t >= t0 && w == 0
                
                t = t0;
                XA = XA-1;
                XB = XB+1;
                w = 1;
                
            else

                % Build a sampling tower

                ir2 = 1;
                r2 = rand;
                cumul(1) = T_A_rep;
                cumul(2) = T_A_rep+T_A_death;
                cumul(3) = T_A_rep+T_A_death+T_B_rep;
                cumul(4) = T;

                % Determine which reaction occurs and compute the number of
                % microbes

                while cumul(ir2) < r2*T

                    ir2 = ir2+1;

                end

                if ir2 == 1

                    XA = XA+1;

                elseif ir2 == 2

                    XA = XA-1;

                elseif ir2 == 3

                    XB = XB+1;

                elseif ir2 == 4

                    XB = XB-1;

                end
            
            end

        end
        
        XB_list(1, i) = XB;
        
    end
    
    % Compute the fixation probability
     
    pfix = length(XB_list(XB_list ~= 0))/Nit;
   
end
