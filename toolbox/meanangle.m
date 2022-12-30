function ANGLEMEAN = meanangle(ANGLES,LIMITS,DIM)

% MEANANGLE  Calculates the mean angle, correcting for angle discontinuities
%   ANGLEMEAN = MEANANGLE(ANGLES,LIMITS,DIM) calculates the mean angle of the
%   input vector/matrix ANGLES taking into account discontinuities. If DIM is
%   specified, the mean over dimension DIM is taken. ANGLES can be up to 2D.
%
%   Synopsis: ANGLEMEAN = MEANANGLE(ANGLES,LIMITS,DIM)
%   Input:    ANGLES   : Input vector or matrix of angles (in degrees). Up
%                        to two dimensions.
%             LIMITS   : Definition of angles: [-180 180] or [0 360] degrees.
%             DIM      : Dimension over which to operate.
%   Output:   ANGLEMEAN: Mean of ANGLES (in degrees).
%
%   Pascal de Theije, v2.0, 13 April 2005
%   Copyright (c) 2005, TNO-D&V
%   All Rights Reserved

%
% Check if input angles match with specified angle limits.
%
if ~isempty(find((ANGLES<min(LIMITS))|(ANGLES>max(LIMITS))))
  disp('meanangle.m: Input angles do not match with limits.')
  return
end

% If the dimension is not specified, operate over the first dimension.
if (nargin == 2)
  DIM = 1;
end

%
% Transpose matrix, if necessary, so that all operations have to be done
% over the second dimension.
%
if (DIM == 1)
  ANGLES = ANGLES.';
end

%
% In order to find the best estimate of the mean angle, calculate the
% in-product of an arbitrary vector [a b] and all specified angles, and find 
% that direction that gives the highest inproduct. The in-product is given by
% [a b].*[sum(cos(ANGLES)) sum(sin(ANGLES))]. The derivative of this w.r.t.
% 'a' is set to zero, which gives the solution a=sqrt(C^2/(C^2+S^2)).
%
C  = sum(cos(ANGLES*pi/180),2);
S  = sum(sin(ANGLES*pi/180),2);
CS = C.^2 + S.^2;
% Calculate vector that gives highest inproduct.
a  = sqrt(C.^2./CS);
b  = sqrt(S.^2./CS);

%
% From the way 'a' and 'b' are solved, they can be either positive or negative,
% while in practice only one of these combinations is the true one. Find the
% combination of 'a' and 'b' that gives the lowest in-product. This angle
% is used to 'split' the angle circle.
%
temp = (C.*a)*[-1 1 -1 1] + (S.*b)*[-1 -1 1 1];
[dummy,ind] = max(temp,[],2);
ind2 = find(ind==1 | ind==3);
a(ind2) = -a(ind2);
ind2 = find(ind==1 | ind==2);
b(ind2) = -b(ind2);

%
% Find angle that should be used to 'split' the angle circle.
%
cut_angle = atan2(b,a) * 180 / pi;

%
% Split the angle circle at 'cut_angle'.
%
ANGLES = ang180(ANGLES - cut_angle*ones(1,size(ANGLES,2)));
% Then take 'normal' mean.
ANGLEMEAN = mean(ANGLES,2);
% Undo splitting angle circle.
ANGLEMEAN = ang180(ANGLEMEAN + cut_angle);

%
% Transform output to angles between 0 and 360 degrees, if necessary.
%
if (LIMITS == [0 360])
  ind = find(ANGLEMEAN < 0);
  ANGLEMEAN(ind) = ANGLEMEAN(ind) + 360;
end

%
% Transpose matrix, if necessary.
%
if (DIM == 1)
  ANGLEMEAN = ANGLEMEAN.';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = ANG180(X);

% ANG180   Unwraps the input angle to an angle between -180 and 180 degrees.
%
%   Pascal de Theije, v1.0, 24 February 1999
%   Copyright (c) 2003, TNO-FEL
%   All Rights Reserved

Y = mod(X-(1e-10)-180,360) - 180 + (1e-10);
                         % The term 1e-10 is to set an input of -180 to 180.
