clear; % clear the workspace
clc; % clear the terminal

%% File Name: Interim Project 1

format compact

%% Variable Initialization

A = 35.3; 
Ayy = 8.555; 
Azz = 23.03; 
Iyy = 495; 
Izz = 1380; 
J = 9.37; 
E = 29000; 
v = 0.3; 
node1 = AJCT_Node([0; 0; 0],1);
node2 = AJCT_Node([80; 80; 0],2);
webdir = [-1, 1, 0]/sqrt(2); 
w = [1;1;1];
element1 = AJCT_Element(A, Ayy, Azz, E, [node1; node2], Izz, Iyy, J, v, webdir, w);

%% Main Code