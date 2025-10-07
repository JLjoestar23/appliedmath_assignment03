% example function 01
% describes first order system drived by sinusoids
function dXdt = rate_func01(t, X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end