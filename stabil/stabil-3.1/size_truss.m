function [S,dSdx] = size_truss(Node,dNodedx)

% SIZE_TRUSS    Compute truss element size (length).
%
%   S = SIZE_TRUSS(Node) computes the element size (length) of a two node
%   truss element.
%
%   [S,dSdx] = SIZE_TRUSS(Node,dNodedx) additionally computes the
%   derivatives of the element size with respect to the design variables x.
%
%   Node        Node definitions               	[x y z] (3 * 3)
%   dNodedx     Node definitions derivatives   	(SIZE(Node) * nVar)
%   S           Element size
%   dSdx        Element size derivatives
%
%   See also ELEMSIZES, ELEMVOLUMES, SIZE_BEAM.

% Wouter Dillen
% December 2017

if nargin<2, dNodedx = []; end

[S,dSdx] = size_beam(Node,dNodedx);


    
    
