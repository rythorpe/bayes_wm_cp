function [LR, UP, PE]=computeLR(outcomes, predictions, newBlock, coordSys)
% Compute update (UP), prediction error (PE), and learning rate (LR=UP/PE)
% for either cartesian or polar coordinates
% Last modified by Ryan Thorpe on 2023.07.03

% Define default coordinate system
if nargin<4
    coordSys = 'cartesian';
end

% find the last trial of each block
a=find(newBlock)-1;
a=a(a>1);

UP=nan(length(outcomes),1);   %nans
%UP_options=nan(length(outcomes),3);

switch coordSys
    case 'cartesian'
        UP(1:end-1) = predictions(2:end)-predictions(1:end-1);
        PE = outcomes-predictions;
    case 'polar'
        PE = circ_dist(outcomes,predictions); % Minimal distance between outcome and prediction orientations
        UP(1:end-1) = PE(1:end-1) + circ_dist(predictions(2:end),predictions(1:end-1)+PE(1:end-1)); % Update always travels in the direction of PE
    case 'polarNoCorrect'
        PE = circ_dist(outcomes,predictions); % Minimal distance between outcome and prediction orientations 
        UP(1:end-1) = circ_dist(predictions(2:end),predictions(1:end-1)); % Distance one way around the circle
    otherwise
        error("computeLR error: didn't recognize coordinate system!! Use coordSys='cartesian' or 'polar'.")
end

% Learing rate
LR=UP./PE;

% Undefined for the first trial of each block
UP(a)=nan;
PE(a)=nan;
LR(a)=nan;
end
