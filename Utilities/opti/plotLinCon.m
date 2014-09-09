function plotLinCon(A,b,Aeq,beq)
%PLOTLINCON Plot Linear Constraints on the current figure
%   plotLinCon(A,b,Aeq,beq)

%   Copyright (C) 2011 Jonathan Currie (I2C2)

xl = xlim; yl = ylim;

hold on;

%Plot Inequality Constraints
if(~isempty(A))
    A = full(A);
    %Check for 1D, fake second column as zeros
    if(size(A,2) == 1)
        A = [A zeros(size(A,1),1)];
    end
    
    for i = 1:length(b)  
        if(any(A(i,:)))
            if(all(A(i,:)) && b(i)) %normal, full A + b
                %%disp('normal')
                c = b(i) / A(i,2);
                m = -(b(i) / A(i,2)) / (b(i) / A(i,1));
                xi = (yl - c) / m;
            elseif(all(A(i,:)) && ~b(i)) %b is zero, line through origin
                %%disp('b is zero, full A')
                c = 0;
                m = -A(i,1)/A(i,2);
                xi = (yl - c) / m;
            elseif(b(i)) %bounds, b scales
                %%disp('partially empty A, b scales')                
                m = inf;
            else %bounds, b is zero
                %%disp('partially empty A, empty b');
                m = inf;
            end

            %Normal Line
            if(~isinf(m) && ~isnan(m))
                %Grab a point (+1,+1) or (-1,+1)
                if(m <= 0)
                    x = [xi(1)+1; yl(1)+1];
                else
                    x = [xi(1)-1; yl(1)+1];
                end
                %Evaluate point to check direction of patch
                in = A(i,:)*x;
                %Patch based on type + gradient
                if(in > b(i)) % <=
                    %%disp('leq')
                    if(m <= 0)
                        %%disp('neg grad')
                        patch([xi(1) xl(2) xl(2) xi(2)],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.3)
                    else
                        %%disp('pos grad')
                        patch([xi(1) xl(1) xl(1) xi(2)],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.3)
                    end
                else % >=
                    %%disp('geq')
                    if(m <= 0)
                        %%disp('neg grad')
                        patch([xi(1) xl(1) xl(1) xi(2)],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.3)
                    else
                        %%disp('pos grad')
                        patch([xi(1) xl(2) xl(2) xi(2)],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.3)
                    end
                end
            %Bounds    
            else
                if(A(i,1) == 0)
                    %%disp('A1 empty')
                    a = A(i,2);
%                     if(xl(1) == -1), xxl = 0.1; else xxl = 1; end
%                     in = a*(xl(1)+xxl) > b(i);
                    i1 = a*(yl(1))-b(i); i2 = a*(yl(2))-b(i);
                    if(i1 < i2), in = false; else in = true; end
                    vbnd = b(i)/a;
                    bndtype = 1; %horizontal (y)
                else
                    %%disp('A2 empty')
                    a = A(i,1);    
%                     if(yl(1) == -1), yyl = 0.1; else yyl = 1; end
%                     in = a*(yl(1)+yyl) > b(i);
                    i1 = a*(xl(1))-b(i); i2 = a*(xl(2))-b(i);
                    if(i1 < i2), in = false; else in = true; end
                    vbnd = b(i)/a;
                    bndtype = 0; %vertical (x)
                end                
                if(bndtype == 0) %vertical
                    %%disp('vertical bound')
                    if(in)    
                        %%disp('lb')
                        patch([vbnd vbnd xl(1) xl(1)],[yl(1) yl(2) yl(2) yl(1)],'y','FaceAlpha',0.3)
                    else
                        %%disp('ub')
                        patch([vbnd vbnd xl(2) xl(2)],[yl(1) yl(2) yl(2) yl(1)],'y','FaceAlpha',0.3)
                    end      
                else %horizontal
                    %%disp('horizontal bound')
                    if(in)        
                        %%disp('lb')
                        patch([xl(1) xl(1) xl(2) xl(2)],[vbnd yl(1) yl(1) vbnd],'y','FaceAlpha',0.3)
                    else
                        %%disp('ub')
                        patch([xl(1) xl(1) xl(2) xl(2)],[vbnd yl(2) yl(2) vbnd],'y','FaceAlpha',0.3)
                    end
                end
            end
        end
    end
end

%Plot Equality Constraints
if(~isempty(Aeq))
    Aeq = full(Aeq);
    %Check for 1D, fake second column as zeros
    if(size(Aeq,2) == 1)
        Aeq = [Aeq zeros(size(Aeq,1),1)];
    end

    for i = 1:length(beq)
        if(any(Aeq(i,:)))
            if(all(Aeq(i,:)) && beq(i)) %normal, full A + b
                %disp('eq normal')
                c = beq(i) / Aeq(i,2);
                m = -(beq(i) / Aeq(i,2)) / (beq(i) / Aeq(i,1));
                xi = (yl - c) / m;
            elseif(all(Aeq(i,:)) && ~beq(i)) %b is zero, line through origin
                %disp('eq b is zero, full A')
                c = 0;
                m = -Aeq(i,1)/Aeq(i,2);
                xi = (yl - c) / m;
            elseif(beq(i)) %bounds, b scales
                %disp('eq partially empty A, b scales')                
                m = inf;
            else %bounds, b is zero
                %disp('eq partially empty A, empty b');
                m = inf;
            end

            %Normal Line
            if(~isinf(m) && ~isnan(m))
                line(xi,yl);                
            %Straight Line (H or V)
            else
                if(Aeq(i,1) == 0)
                    %disp('eq A1 empty')
                    a = Aeq(i,2);
                    vbnd = beq(i)/a;
                    line(xl,[vbnd vbnd]);
                else
                    %disp('eq A2 empty')
                    a = Aeq(i,1);                    
                    vbnd = beq(i)/a;
                    line([vbnd vbnd],yl);
                end

            end
        end
    end
end

hold off;

end

