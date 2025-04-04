% ������Դ��https://www.researchgate.net/publication/371247711_A_mixed_strategy_improved_dung_beetle_optimization_algorithm_and_its_application
%%% 3���Ľ����ԣ�1.�ѵ㼯��ʼ������ 2.�������������� 3.Levy���в���

function [fMin , bestX, Convergence_curve ] = MSDBO(N, Max_iteration,lb,ub,dim,fobj)
        
P_percent = 0.2;    % The population size of producers accounts for "P_percent" percent of the total population size       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pNum = round( N *  P_percent );    % The population size of the producers   


lb= lb.*ones( 1,dim );    % Lower limit/bounds/     a vector
ub= ub.*ones( 1,dim );    % Upper limit/bounds/     a vector



%Initialization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% �Ľ�����һ���ѵ㼯��ʼ������ %%%
GD=initialization_gd(N,dim,lb,ub);
% % ����GD��Ϊ���ļѵ㼯
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : N
    
    x( i, : ) = GD(i,:);  
    fit( i ) = fobj( x( i, : ) ) ;                       
end

pFit = fit;                       
pX = x; 
 XX=pX;    
[ fMin, bestI ] = min( fit );      % fMin denotes the global optimum fitness value
bestX = x( bestI, : );             % bestX denotes the global optimum position corresponding to fMin

 % Start updating the solutions.
for t = 1 : Max_iteration    
       
        [fmax,B]=max(fit);
        worse= x(B,:);   
       r2=rand(1);
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : pNum    
        if(r2<0.9)
            r1=rand(1);
          a=rand(1,1);
          if (a>0.1)
           a=1;
          else
           a=-1;
          end
          b = rand();
    x( i , : ) =  pX(  i , :)+0.3*abs(pX(i , : )-worse)+a*0.1*(XX( i , :)); % Equation (1)
       else
            
           aaa= randperm(180,1);
           if ( aaa==0 ||aaa==90 ||aaa==180 )
            x(  i , : ) = pX(  i , :);   
           end
         theta= aaa*pi/180;   
       
       x(  i , : ) = pX(  i , :)+tan(theta).*abs(pX(i , : )-XX( i , :));    % Equation (2)      

        end
      
        x(  i , : ) = Bounds( x(i , : ), lb, ub );    
        fit(  i  ) = fobj( x(i , : ) );
    end 
 [ fMMin, bestII ] = min( fit );      % fMin denotes the current optimum fitness value
  bestXX = x( bestII, : );             % bestXX denotes the current optimum position 

 R=1-t/Max_iteration;                           %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Xnew1 = bestXX.*(1-R);
 Xnew2 =bestXX.*(1+R);                    %%% Equation (3)
 Xnew1= Bounds( Xnew1, lb, ub );
 Xnew2 = Bounds( Xnew2, lb, ub );
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Xnew11 = bestX.*(1-R);
 Xnew22 =bestX.*(1+R);                     %%% Equation (5)
 Xnew11= Bounds( Xnew11, lb, ub );
 Xnew22 = Bounds( Xnew22, lb, ub );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% �Ľ����Զ����������������� %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 5;
l = 2*rand()-1;
z = exp(k*cos(pi*(t/Max_iteration)));

for i = ( pNum + 1 ) :12                  % Equation (4)
    % x( i, : )=bestXX+((rand(1,dim)).*(pX( i , : )-Xnew1)+(rand(1,dim)).*(pX( i , : )-Xnew2));
    x( i, : )= bestXX+(exp(z*l) * cos(2*pi*l) * (rand(1,dim)).*(pX( i , : )-Xnew1)+(exp(z*l) * cos(2*pi*l) * rand(1,dim)).*(pX( i , : )-Xnew2));  % �˴��иĽ�������ɱ���������
    x(i, : ) = Bounds( x(i, : ), Xnew1, Xnew2 );
    fit(i ) = fobj(  x(i,:) ) ;
end
   
   for i = 13: 19                  % Equation (6)


       % x( i, : )=pX( i , : )+((randn(1)).*(pX( i , : )-Xnew11)+((rand(1,dim)).*(pX( i , : )-Xnew22)));
       x( i, : ) =  exp(z*l) * cos(2*pi*l) * pX( i , : )+((randn(1)).*(pX( i , : )-Xnew11)+((rand(1,dim)).*(pX( i , : )-Xnew22)));  %�иĽ���������
       x(i, : ) = Bounds( x(i, : ),lb, ub);
       fit(i ) = fobj(  x(i,:) ) ;

   end
  
  for j = 20 : N                 % Equation (7)
       x( j,: )=bestX+randn(1,dim).*((abs(( pX(j,:  )-bestXX)))+(abs(( pX(j,:  )-bestX)))).*0.5;
      x(j, : ) = Bounds( x(j, : ), lb, ub );
      fit(j ) = fobj(  x(j,:) ) ;
  end
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% �Ľ���������Levy���� %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Temp=bestX+1+Levy(dim); % �ο������е�(13)ʽ�����Ե�ǰ���Ž�ִ��
  Temp(Temp>ub) = ub(Temp>ub);
  Temp(Temp<lb) = lb(Temp<lb);
  fitvalue = fobj(Temp);
  if(fitvalue <fit(i))
      x( i, : ) = Temp;
      fit(i) = fitvalue;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  % Update the individual's best fitness vlaue and the global best fitness value
     XX=pX;
    for i = 1 : N 
        if ( fit( i ) < pFit( i ) )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end
        
        if( pFit( i ) < fMin )
           % fMin= pFit( i );
           fMin= pFit( i );
            bestX = pX( i, : );
          %  a(i)=fMin;
            
        end
    end
  
     Convergence_curve(t)=fMin;
  
    
  
end
end

% Application of simple limits/bounds
function s = Bounds( s, Lb, Ub)
  % Apply the lower bound vector
  temp = s;
  I = temp < Lb;
  temp(I) = Lb(I);
  
  % Apply the upper bound vector 
  J = temp > Ub;
  temp(J) = Ub(J);
  % Update this new move 
  s = temp;
%---------------------------------------------------------------------------------------------------------------------------
end

%% ����Levy������������Ӻ���
function o=Levy(D)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,D)*sigma;
v=randn(1,D);
step=u./abs(v).^(1/beta);
o=0.01*step;
end

% This function initialize the first population of search agents
function Positions=initialization_gd(SearchAgents_no,dim,lb,ub)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions = Goodnode(SearchAgents_no,dim).*(ub-lb)+lb;
    Positions = Bounds(Positions,lb,ub);   % �߽紦��
end

% If each variable has a different lb and ub
if Boundary_no>1
    GoodnodeValue = Goodnode(SearchAgents_no,dim);
    Positions = GoodnodeValue.*(ub-lb)+lb;
    Positions = Bounds(Positions,lb,ub);   % �߽紦��
end
end

function [GD] = Goodnode(num,dim)
% num�ǵ���
% dim��ά��
tmp1 = [1:num]'*ones(1, dim);
Ind = [1:dim];
prime1 = primes(100*dim);
[p,q]=find(prime1 >= (2*dim+3));
tmp2 = (2*pi.*Ind)/prime1(1,q(1));
tmp2 = 2*cos(tmp2);
tmp2 = ones(num,1)*tmp2;
GD = tmp1.*tmp2;
GD = mod(GD,1);
end
