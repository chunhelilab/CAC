
function [out, actions, occupancies]=IM_OLACFunctional(actions, FP)
 
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
Epsilon = 200000; 
rates=exp(-actions/Epsilon);
rates(logical(eye(sz))) = -sum(rates,2); 
occupancies = null(rates')/sum(null(rates'));
%This calculates the occupancies in a canonical way by finding the steady
%state of the continuous time Markov chain. 
%%%%%5 stable states as we maximize the occupancy of A state,4 is the A state.
 s=zeros(5,1);
 s(4)=1;
 s(1)=s(4)*(rates(1,4)-rates(4,1));
 s(2)=s(4)*(rates(2,4)-rates(4,2));
 s(3)=s(4)*(rates(3,4)-rates(4,3));
 out=-sum(s)

end