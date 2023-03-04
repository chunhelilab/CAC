%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MISE_OLACFunctional.
%This is a function to evaluate a matrix of actions, along with the corresponding
%fixed point. 
%Here we evaluate the occupancy of the third state of the system. 
%Note that since we want to maximize this function, we actually minimize
%its negative, as is standard practice. 
%This function can account for when other states have disappeared. 


%Author: Daniel K. Wells (? 2015
%Ver 1.0
%Email: dannykwells@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [out, actions, occupancies]=IM_OLACFunctional(actions, FP)
%For our test functional, we optimize the occupancy of the third fixed
%point in the MISE system, assuming a noise level of epsilon = 0.02;


%Find those transitions which are faster when indirect. 
sz = size(actions, 2);
for i=1:sz
    for j=1:sz
        for k=1:sz
            if actions(i,j) + actions(j,k) < actions(i,k)
                actions(i,k) = Inf; 
            end
        end
    end
end

actions(logical(eye(sz))) = Inf; 
%Epsilon = 200000; %CL
Epsilon =20000;
rates=exp(-actions/Epsilon);
rates(logical(eye(sz))) = -sum(rates,2); 
occupancies = null(rates')/sum(null(rates'))
%This calculates the occupancies in a canonical way by finding the steady
%state of the continuous time Markov chain. 

%This is the occupancy of the third state. 
 res = -sum(occupancies(FP==4)) 

%out = -sum(occupancies(FP==2)); %%CL, occupancy of the second state, Immune state 

%%%%%%%%%calculate A-M states
%{
 s=zeros(4,1);
 s(4)=1;
 s(1)=s(4)*(rates(1,4)-rates(4,1));
 s(2)=s(4)*(rates(2,4)-rates(4,2));
 s(3)=s(4)*(rates(3,4)-rates(4,3));
 out=-sum(s)
%}
%%%%%%%%%calculate A states
%
 s=zeros(5,1);
 s(5)=1;
 s(1)=s(5)*(rates(1,5)-rates(5,1));
 s(2)=s(5)*(rates(2,5)-rates(5,2));
 s(3)=s(5)*(rates(3,5)-rates(5,3)-rates(3,4));
 s(4)=s(5)*(rates(4,5)-rates(5,4));
 out=-sum(s(1:4))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

end