% Code created by Loïc Marrec

function tauaf = Gillespie_fct(Nit, n, theta, gA, XA_i, gB, XB_i, K, mu)

    t_list = NaN(Nit, 1);
    Napp = NaN(Nit, 1);
    q = 1;

    while q <= Nit

        % Initialization

        XA = XA_i;              % Initialization of the number of A individuals
        XB = XB_i;              % Initialization of the number of B individuals
        XBlist = nan(1, 5e1);
        tapplist = nan(1, 5e1);
        t = 0;                 % Initialization of time
        cumul = zeros(5, 1);    % To build the sampling tower
        qapp = 1;

        while XA ~= 0 
            
            % Compute the wild-type fitness and the mutant fitness

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
                    
                tapplist(1, qapp) = t;
                XBlist(1, qapp) = 1;
                qapp = qapp+1;

            elseif ir2 == 3

                XA = XA-1;

            elseif ir2 == 4

                XB = XB+1;
                    
                ir3 = 1;
                r3 = rand;
                cumul2 = cumsum(XBlist(1 : length(XBlist)-length(XBlist(isnan(XBlist)))));
                    
                while cumul2(ir3) < r3*cumul2(length(cumul2))

                     ir3 = ir3+1;

                end
                    
                XBlist(1, ir3) = XBlist(1, ir3)+1;
                    

            elseif ir2 == 5

                XB = XB-1;
                    
                ir4 = 1;
                r4 = rand;
                cumul3 = cumsum(XBlist(1 : length(XBlist)-length(XBlist(isnan(XBlist)))));
                    
                while cumul3(ir4) < r4*cumul3(length(cumul3))

                    ir4 = ir4+1;

                end
                    
                XBlist(1, ir4) = XBlist(1, ir4)-1;

            end
            
        end
        
        if XA == 0 && XB ~= 0
            
            [~, ind] = max(XBlist);
            
            t_list(q) = tapplist(ind);
            Napp(q) = qapp;
            q = q+1;
            
        end
        
    end
    
    tauaf = nanmean(t_list);
   
end
