% Dung Beetle Optimizer: (DBO) (demo)
% Programmed by Jian-kai Xue    
% Updated 28 Nov.2022. 

function [fMin, bestX, Convergence_curve] = DBO(pop, M, c, d, dim, fobj)

P_percent = 0.2; % Proportion of dung beetles performing ball-rolling behavior   
pNum = round(pop * P_percent); % Number of dung beetles performing ball-rolling behavior  
lb = c.*ones(1,dim);    % lower bound vector for the value
ub = d.*ones(1,dim);    % upper bound vector for the value

% population initialization
for i = 1 : pop
    x(i,:) = lb + (ub - lb) .* rand(1,dim);  
    fit(i) = fobj(x(i,:)) ;                       
end
pFit = fit;                       
pX = x; 
XX = pX;    
[fMin, bestI] = min(fit); % fMin is the value of global optimal fitness
bestX = x(bestI,:);       % bestX is the global optimal solution corresponding to fMin

% Beginning iteration for individual updates
for t = 1 : M 

[fmax, B] = max(fit);
worse = x(B,:);   
r2 = rand(1);

% Rolling Ball Behavior Dung Beetle Position Updates (divided into No Obstacle Mode and Obstacle Mode)
for i = 1 : pNum    
    if(r2<0.9)
        % ① No Obstacle Mode
        a = rand(1,1);
        if (a>0.1)
            a = 1;
        else
            a = -1;
        end
        x(i,:) = pX(i,:)+0.3*abs(pX(i,:)-worse)+a*0.1*(XX(i,:));
    else
        % ② Obstacle Mode
        aaa = randperm(180,1);
        if (aaa==0 ||aaa==90 ||aaa==180)
            x(i,:) = pX(i,:);   
        end
        theta = aaa*pi/180;   
        x(i,:) = pX(i,:)+tan(theta).*abs(pX(i,:)-XX(i,:)); 
    end
    % 越界校正
    x(i,:) = Bounds(x(i,:), lb, ub);    
    fit(i) = fobj(x(i,:));
end

[fMMin, bestII] = min(fit);   % fMin is the current optimal fitness value
bestXX = x(bestII,:);         % bestXX is the current optimum solution
R = 1-t/M;                         

Xnew1 = bestXX.*(1-R); 
Xnew2 = bestXX.*(1+R);                    
Xnew1 = Bounds(Xnew1, lb, ub);
Xnew2 = Bounds(Xnew2, lb, ub);

% Reproductive Behavior of Dung Beetles Location Updates
for i = (pNum + 1) : 12               
    x(i,:) = bestXX+((rand(1,dim)).*(pX(i,:)-Xnew1)+(rand(1,dim)).*(pX(i,:)-Xnew2));
    x(i,:) = Bounds(x(i,:), Xnew1, Xnew2);
    fit(i) = fobj(x(i,:)) ;
end

Xnew11 = bestX.*(1-R); 
Xnew22 = bestX.*(1+R);                    
Xnew11 = Bounds(Xnew11, lb, ub);
Xnew22 = Bounds(Xnew22, lb, ub);

% Foraging Behavior of Small Dung Beetles Location Updates
for i = 13 : 19        
    x(i,:) = pX(i,:)+((randn(1)).*(pX(i,:)-Xnew11)+((rand(1,dim)).*(pX(i,:)-Xnew22)));
    x(i,:) = Bounds(x(i,:), lb, ub);
    fit(i) = fobj(x(i,:));
end

% Stealing Behavior Dung Beetle Location Updates
for j = 20 : pop       
    x(j,:) = bestX+randn(1,dim).*((abs((pX(j,:)-bestXX)))+(abs((pX(j,:)-bestX))))./2;
    x(j,:) = Bounds(x(j,:), lb, ub);
    fit(j) = fobj(x(j,:)) ;
end

% Update the individual optimum and the global optimum
XX = pX;
for i = 1 : pop 
    if (fit(i) < pFit(i))
        pFit(i) = fit(i);
        pX(i,:) = x(i,:);
    end

    if( pFit(i) < fMin)
        fMin = pFit(i);
        bestX = pX(i,:);
    end
end

% Iteration Curve Records
Convergence_curve(t) = fMin;
end

% Out of bounds correction function
function s = Bounds(s, Lb, Ub)
temp = s;
I = temp < Lb;
temp(I) = Lb(I);

J = temp > Ub;
temp(J) = Ub(J);
s = temp;

