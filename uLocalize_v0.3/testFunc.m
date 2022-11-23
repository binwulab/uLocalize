function [O1, O2, O3]=testFunc(I1,I2,I3)
disp(['# of Input: ', num2str(nargin)]);
O1=I1;
O2=I2;
O3=I3;
disp(['# of Output: ', num2str(nargout)]);
end