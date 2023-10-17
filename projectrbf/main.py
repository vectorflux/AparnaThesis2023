# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np
#import cupy as cp
#from fixtures.context import *
#import ghex.unstructured as unstructured

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    myzero = np.zeros([2,3])
    print(myzero)
    length = np.shape(myzero)
    print(length[0], length[1])



# See PyCharm help at https://www.jetbrains.com/help/pycharm/

#
# Inputs:
# % d = dimension of state space (positive integer)
# % k = smoothness parameter (non-negative integer)
# %
# % Outputs:
# % out = function handle to the Wendland RBF
# %       accepts values in [0,inf)
# %       normalised such that out(0) = 1
# %
# % Background:
# % Wendland(d,k) is a polynomial on [0,1] and vanishes on (1,inf)
# % Wendland(d,k) is C^2k(R); i.e. it has 2k continuous derivatives
# % As a RBF, Wendland(d,k) is positive definite on R^d
# %
# % Example:
# % Wendland(d,0) = (1-r)_+^l
# % Wendland(d,1) = (1-r)_+^(l+1) * ((l+1)r + 1)
# % Wendland(d,2) = (1-r)_+^(l+2) * ((l^2+4l+3)r^2 + (3l+6)r + 3)
# %                 in each case up to a multiplicative constant
# %                 z_+^l := max(0,z)^l
# %                 l = floor(d/2) + k + 1
# %
# % © Chris Oates 2017.function out = Wendland(d,k)
# l = floor(d/2) + k + 1;
# phi = @(r) (1-r)^l; % power function (Sec 11.2 in Fasshauer)
# for i = 1:k
#     phi = I(phi); % iterative integration (Defn 11.2 in Fasshauer)
# end
# nor = phi(0);
# syms r
# out = phi(r) / nor; % normalise at r = 0
# out = heaviside(r) * heaviside(1-r) * out; % truncation to compact support
# out = matlabFunction(out);
# end% Integration operator
# % f(r) -> \int_r^\infty t f(t) dt
# function out = I(in)
# syms r
# integrand = @(t) t * in(t);
# out = int(integrand,r,1);
# out = matlabFunction(out);
# end
