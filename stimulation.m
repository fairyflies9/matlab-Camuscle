function [stim] = stimulation(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to determine if stimulus is on (1) or off (0).
%Inputs
% t = time in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(t<0.3),
    stim=1;
else
    stim=0;
end

end
