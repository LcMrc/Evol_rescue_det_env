% Code created by Loïc Marrec

function pr = Gillespie_fct(Nit, n, theta, gA, XA_i, gB, XB_i, K, mu)

    XB_list = NaN(1, Nit);

    for i = 1 : Nit

        % Initialization

        XA = XA_i;              % Initialization of the number of A individuals
        XB = XB_i;              % Initialization of the number of B individuals
        t = 0;                  % Initialization of time
        cumul = zeros(5, 1);    % To build the sampling tower

        while XA ~= 0 

            fA = sigm(t, theta, n);
            fB = 1-sigm(t, theta, n);
            
            % Compute the transition rates

            T_A_rep_no_mut = fA*(1-(XA+XB)/K)*(1-mu)*XA; 
            T_A_rep_mut = fA*(1-(XA+XB)/K)*mu*XA;
            T_A_death = gA*XA;
            T_B_rep = fB*(1-(XA+XB)/K)*XB;
            T_B_death = gB*XB;
            T = T_A_rep_no_mut+T_A_rep_mut+T_A_death+T_B_rep+T_B_death;

            % Compute tau and then the new time

            r1 = rand;
            tau = 1/T*log(1/r1);
            t = t+tau;

                % Build a sampling tower

                ir2 = 1;
                r2 = rand;
                cumul(1) = T_A_rep_no_mut;
                cumul(2) = T_A_rep_no_mut+T_A_rep_mut;
                cumul(3) = T_A_rep_no_mut+T_A_rep_mut+T_A_death;
                cumul(4) = T_A_rep_no_mut+T_A_rep_mut+T_A_death+T_B_rep;
                cumul(5) = T;

                % Determine which reaction occurs and compute the number of
                % individuals

                while cumul(ir2) < r2*T

                    ir2 = ir2+1;

                end

                if ir2 == 1

                    XA = XA+1;
                    
                elseif ir2 == 2
                    
                    XB = XB+1;

                elseif ir2 == 3

                    XA = XA-1;

                elseif ir2 == 4

                    XB = XB+1;

                elseif ir2 == 5

                    XB = XB-1;

                end
            
        end

        XB_list(1, i) = XB;
        
    end
    
    % Compute the rescue probability
     
    pr = length(XB_list(XB_list ~= 0))/Nit;
   
end
