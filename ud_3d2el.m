function [DEFL,REACT,ELE_FOR,AFLAG,APRATIOS,LIMIT_STATE] = ud_3d2el(...
		nnodes,coord,concen,fixity,nele,ends,A,Izz,Iyy,J,Cw,IsSym,Ysc,Zsc,Betay,Betaz,Betaw,Zzz,Zyy,Ayy,Azz,...
		E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,E_t,anatype,sol_scheme,numsteps,...
		ratio_req,stop_ratio,restart,defl,react,ele_for,...
		apratios,limit_state,h_stat_mes)
%UD_3D2EL performs a user defined three-dimensional
% second-order elastic analysis of a structural system.
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Functions Called
%              < to be defined by the student >
%
%  Dictionary of Variables
%     Input Information:
%       nnodes         ==  total number of nodes
%       coord(i,1:3)   ==  node i's coordinates
%                            coord(i,1) = X coordinate
%                            coord(i,2) = Y coordinate
%                            coord(i,3) = Z coordinate
%       concen(i,1:6)  ==  concentrated loads for node i's 6 d.o.f.
%                            concen(i,1) = force in global X direction
%                            concen(i,2) = force in global Y direction
%                            concen(i,3) = force in global Z direction
%                            concen(i,4) = moment about global X axis
%                            concen(i,5) = moment about global Y axis
%                            concen(i,6) = moment about global Z axis
%       fixity(i,1:6)  ==  prescribed displacements for node i's 6 d.o.f.
%                          Note: A free d.o.f. will have a value of NaN
%                          Examples: If fixity(15,3) = NaN, then node 15's
%                                      Z-disp component is free;
%                                    If fixity(2,6) = 0.0, then node 2's
%                                      Z-rotation component is supported;
%                                    If fixity(5,2) = -2.1, then node 5's
%                                      Y-disp component is supported and defined
%                                      with a settlement of -2.1 units.
%                            fixity(i,1) = prescribed disp. in global X direction
%                            fixity(i,2) = prescribed disp. in global Y direction
%                            fixity(i,3) = prescribed disp. in global Z direction
%                            fixity(i,4) = prescribed rotation about global X axis
%                            fixity(i,5) = prescribed rotation about global Y axis
%                            fixity(i,6) = prescribed rotation about global Z axis
%       nele           ==  total number of elements
%       ends(i,1:14)   ==  element i's nodal information
%                            ends(i,1) = start node #
%                            ends(i,2) = finish node #
%                            ends(i,3) = flag to indicate whether or not flexural
%                            moments are released at start node.  ends(i,3)=0 both not
%                            released (rigid connection); ends(i,3)=1 both flexural
%                            moments are released (pinned connection); ends(i,3)=2
%                            at least one of the flexural moments are partially or fully
%                            released (see below for connection stiffness attributes)
%                            ends(i,4) = flag to indicate whether or not flexural
%                            moments are released at finish node.  ends(i,4)=0 both not
%                            released (rigid connection); ends(i,4)=1 both flexural
%                            moments are released (pinned connection); ends(i,4)=2
%                            at least one of the flexural moments are partially or fully
%                            released (see below for connection stiffness attributes)
%                            ends(i,5) = flag to indicate the degree of warping
%                            restraint at start node.  ends(i,5)=0 warping free;
%                            ends(i,5)=1 warping fixed; ends(i,5)=2 warping continuous
%                            ends(i,6) = flag to indicate the degree of warping
%                            restraint at finish node.  ends(i,6)=0 warping free;
%                            ends(i,6)=1 warping fixed; ends(i,6)=2 warping continuous
%                            ends(i,7) = rotational spring stiffness at the start
%                            node and about element i's local z-z axis.
%                            ends(i,8) = rotational spring stiffness at the start
%                            node and about element i's local y-y axis.
%                            ends(i,9) = rotational spring stiffness at the finish
%                            node and about element i's local z-z axis.
%                            ends(i,10) = rotational spring stiffness at the finish
%                            node and about element i's local y-y axis.
%                            ends(i,11) = connection moment capacity Mpz at the start
%                            node and about element i's local z-z axis.
%                            ends(i,12) = connection moment capacity Mpy at the start
%                            node and about element i's local y-y axis.
%                            ends(i,13) = connection moment capacity Mpz at the finish
%                            node and about element i's local z-z axis.
%                            ends(i,14) = connection moment capacity Mpy at the finish
%                            node and about element i's local y-y axis.
%       A(i)           ==  element i's cross sectional area
%       Izz(i)         ==  element i's moment of inertia about its local z-z axis
%       Iyy(i)         ==  element i's moment of inertia about its local y-y axis
%       J(i)           ==  element i's torsional constant
%       Cw(i)          ==  element i's warping constant
%       Zzz(i)         ==  element i's plastic section modulus about its local z-z axis
%       Zyy(i)         ==  element i's plastic section modulus about its local y-y axis
%       Ayy(i)         ==  element i's effective shear area along its local y-y axis
%       Azz(i)         ==  element i's effective shear area along its local z-z axis
%       E(i)           ==  element i's material elastic modulus, Young's Modulus
%       v(i)           ==  element i's material Poisson's ratio
%       Fy(i)          ==  element i's material yield strength
%       YldSurf(i)     ==  element i's yield surface maximum values
%                              YldSurf(i,1) = maximum P/Py value
%                              YldSurf(i,2) = maximum Mz/Mpz value
%                              YldSurf(i,3) = maximum My/Mpy value
%       Wt(i)          ==  element i's material weight density
%       webdir(i,1:3)  ==  element i's unit web vector.  This is a unit vector
%                          that defines the element's local y-y axis with respect
%                          to the global coordinate system.  It is based on the
%                          structure's undeformed geometry.
%                              webdir(i,1) = x component of element's unit web vector
%                              webdir(i,2) = y component of element's unit web vector
%                              webdir(i,3) = z component of element's unit web vector
%                          NOTE: An element's 3x3 rotation matrix, [g], is constructed
%                          as follows: First, calculate a unit vector, x_vect, that
%                          describes the element's local x-axis. Second, take the
%                          cross product of x_vect and webdir(i,:) to obtain z_vect,
%                          i.e. z_vect = cross(x_vect,webdir(i,:)). Third, set z_vect 
%                          to a unit vector, i.e. z_vect = z_vect/norm(z_vect).
%                          Finally, the first row of [g] is x_vect, its second row is
%                          webdir(i,:), and its third row is z_vect.
%       beta_ang(i)    ==  element i's web rotation angle.  These values are
%                          provided for those students who are required to calculate
%                          their own unit web vectors (see above).  It is based
%                          on the structure's undeformed geometry.
%                          Note:  MASTAN2 uses the following convention for
%                                 defining a member's default web orientation:
%                                 A vector defing the element's local y-axis
%                                 with respect to the global coordinate system
%                                 will have a positive component in the global
%                                 Y direction.  If the element's local x-axis,
%                                 its length axis, is aligned with the global Y
%                                 axis, then element's local y-axis is aligned
%                                 with global negative X axis.  After this initial
%                                 orientation, element i may be rotated about
%                                 its local x-axis by the amount defined by
%                                 its web rotation angle, beta_ang(i).  The
%                                 angle is in radians and assumes a right-hand
%                                 convention about the local x-axis which runs from
%                                 the element's start node to its finish node.
%       w(i,1:3)         ==  element i's uniform load which references its
%                            local coordinate system
%                              w(i,1) = x component of uniform load
%                              w(i,2) = y component of uniform load
%                              w(i,3) = z component of uniform load
%       thermal(i,1:4)   ==  element i's thermal strain effects which reference its
%                            local coordinate system
%                              thermal(i,1) = coefficient of thermal expansion
%                              thermal(i,2) = change in temperature at centroid
%                              thermal(i,3) = linear temperature gradient in local y-dir
%                                           = (T_up_y - T_btm_y) / depth_y
%                              thermal(i,4) = linear temperature gradient in local z-dir
%                                           = (T_up_z - T_btm_z) / width_z
%       truss            ==  flag to indicate if a truss or not
%                              truss = 0   System is not a truss
%                              truss = 1   System is a truss
%       anatype          ==  flag to indicate which type of analysis is requested
%                              anatype = 1  First-Order Elastic
%                              anatype = 2  Second-Order Elastic
%                              anatype = 3  First-Order Inelastic
%                              anatype = 4  Second-Order Inelastic
%                              anatype = 5  Elastic Buckling (Eigenvalue)
%                              anatype = 6  Inelastic Buckling (Eigenvalue)
%       sol_scheme       ==  flag to indicate which type ofsolution scheme is requested
%                              sol_scheme = 1   Simple Step (Euler)
%                              sol_scheme = 2   Predictor-Corrector (Midpoint R-K)
%                              sol_scheme = 3   Other 1
%                              sol_scheme = 4   Other 2
%                              sol_scheme = 5   Other 3
%       numsteps         ==  requested maximum number of load steps or increments
%       ratio_req        ==  requested load step or increment size
%       stop_ratio       ==  requested maximum applied load ratio
%       restart          ==  flag to indicate if a new or continuing analysis
%                              restart = 1  Start a new analysis
%                              restart = 2  Continue with previous analysis results
%       defl(i,1:6,n)    ==  node i's calculated 6 d.o.f. deflections for all
%                            n previous steps.  See DEFL for additional details.
%       react(i,1:6,n)   ==  node i's calculated 6 d.o.f. reactions for all
%                            n previous steps. See REACT for additional details.
%       ele_for(i,1:12,n)==  element i's calculated 12 element end forces for all
%                            n previous steps.  See ELE_FOR for additional details.
%       apratios(i,1)    ==  applied load ratio at previous step i
%                            Note that apratios is a COLUMN vector.
%                            size(apratios,1) = total number of previous load steps.
%                            See APRATIOS for additional details.
%       limit_state      ==  flag indicating post-limit state analysis
%                              limit_state = 0  pre-limit state, system is loading
%                              limit_state = 1  post-limit state, system is unloading
%       h_stat_mes       ==  handle to graphics object which displays status message.
%                            Example:  % A message at the end of each load step is
%                                      % often informative to the user.  The following
%                                      % Matlab code displays a message which
%                                      % includes the contents of variable istep
%                                      text_mess = ['Performing step #',num2str(istep)];
%                                      set(h_stat_mes,'String',text_mess); drawnow;
%
%     Local Information:
%              < to be defined by the student >
%
%     Output Information:
%       DEFL(i,1:6,m)    ==  node i's calculated 6 d.o.f. deflections for all m steps
%                            of the analysis.  If this is a new analysis, DEFL should be
%                            set equal to [].  If this is a continuing analysis, DEFL
%                            should be first set equal to defl, i.e. DEFL = defl;
%                            For a supported d.o.f., the proper percent of the prescribed
%                            value, fixity(i,dof), should be inserted in DEFL.
%                            Example: DEFL(12,3,5) is node 12's Z-disp. component
%                                     at the end of the fifth load increment.
%                                     If the Z-disp. d.o.f. for node 12
%                                     is supported, then the appropriate percent of
%                                     the prescribed value, fixity(12,3), should
%                                     inserted at DEFL(12,3,5).
%                              DEFL(i,1,m) = displacement in X direction at end of step m
%                              DEFL(i,2,m) = displacement in Y direction at end of step m
%                              DEFL(i,3,m) = displacement in Z direction at end of step m
%                              DEFL(i,4,m) = rotation about X direction at end of step m
%                              DEFL(i,5,m) = rotation about Y direction at end of step m
%                              DEFL(i,6,m) = rotation about Z direction at end of step m
%       REACT(i,1:6,m)   ==  node i's calculated 6 d.o.f. reactions for all m steps
%                            of the analysis.  If this is a new analysis, REACT should be
%                            set equal to [].  If this is a continuing analysis, REACT
%                            should be first set equal to react, i.e. REACT = react;
%                            For free d.o.f., a value of 0.0 should be inserted in REACT. 
%                            Example: REACT(5,6,8) is node 5's Z-moment reaction
%                                     component at the end of the eighth load
%                                     increment.  If the Z-rotation d.o.f. for
%                                     node 5 is free, then a value of 0.0 should
%                                     inserted at REACT(5,6,8).
%                              REACT(i,1,m) = force in X direction at end of step m
%                              REACT(i,2,m) = force in Y direction at end of step m
%                              REACT(i,3,m) = force in Z direction at end of step m
%                              REACT(i,4,m) = moment about X direction at end of step m
%                              REACT(i,5,m) = moment about Y direction at end of step m
%                              REACT(i,6,m) = moment about Z direction at end of step m
%       ELE_FOR(i,1:1?,m)==  element i's calculated element end forces for all m steps
%                            of the analysis.  If this is a new analysis, ELE_FOR should
%                            be set equal to [].  If this is a continuing analysis, 
%                            ELE_FOR should be first set equal to ele_for,
%                            i.e. ELE_FOR = ele_for;
%                            Note: All values reference the element's local
%                                  coordinate system.
%                            Example: ELE_FOR(17,9,5) is element 17's Z-shear
%                                     force at the end node of the element at the
%                                     end of the fifth load increment.
%                              ELE_FOR(i,1,m)  = x-force at start node at end of step m
%                              ELE_FOR(i,2,m)  = y-force at start node at end of step m
%                              ELE_FOR(i,3,m)  = z-force at start node at end of step m
%                              ELE_FOR(i,4,m)  = x-moment at start node at end of step m
%                              ELE_FOR(i,5,m)  = y-moment at start node at end of step m
%                              ELE_FOR(i,6,m)  = z-moment at start node at end of step m
%                              ELE_FOR(i,7,m)  = x-force at end node at end of step m
%                              ELE_FOR(i,8,m)  = y-force at end node at end of step m
%                              ELE_FOR(i,9,m)  = z-force at end node at end of step m
%                              ELE_FOR(i,10,m) = x-moment at end node at end of step m
%                              ELE_FOR(i,11,m) = y-moment at end node at end of step m
%                              ELE_FOR(i,12,m) = z-moment at end node at end of step m
%                            If you are not programming warping torsion, the ELE_FOR
%                            array needs to contain only 12 columns, i.e. ELE_FOR(i,1:12,m)                           
%                            For those programming warping torsion, the bimoments and
%                            rates of twist should be stored as follows.
%                              ELE_FOR(i,13,m) = bimoment at start node at end of step m
%                              ELE_FOR(i,14,m) = bimoment at end node at end of step m
%                              ELE_FOR(i,15,m) = rate of twist at start node at end of step m
%                              ELE_FOR(i,16,m) = rate of twist at end node at end of step m
%       AFLAG            ==  flag to indicate if a successful analysis has or has not been
%                            completed
%                              AFLAG = 1     Successful
%                              AFLAG = 0     Unstable Structure
%                              AFLAG = -1    Analysis Halted: Limit Load Reached
%                              AFLAG = 99    Analysis Halted: No Convergence
%                              AFLAG = inf   No analysis code available
%       APRATIOS(i,1)    ==  applied load ratio at end of load step i
%                            Note that APRATIOS is a COLUMN vector.
%                            size(APRATIOS,1) = total number of load steps completed.
%                            If this is a new analysis, APRATIOS should be
%                            set equal to [].  If this is a continuing analysis, APRATIOS
%                            should be first set equal to apratios,
%                            i.e. APRATIOS = apratios;
%       LIMIT_STATE      ==  flag indicating post-limit state analysis
%                              LIMIT_STATE = 0  limit state not reached, system is loading
%                              LIMIT_STATE = 1  limit state reached or exceeded, system is
%                                               unloading
%
%       Version 1.0/Student's Initials/Date of Modification
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Start by defining all output arrays
%
	if restart == 1  |  isempty(apratios)
		DEFL=[]; REACT=[]; ELE_FOR=[]; APRATIOS=[]; LIMIT_STATE=0;
	else
		DEFL = defl;
		REACT = react;
		ELE_FOR = ele_for;
		APRATIOS = apratios;
		LIMIT_STATE = limit_state;
	end
%
	AFLAG = inf;
%
%  STUDENT NOTE:
%     In order for this routine to become fully active AFLAG
%     must be changed.
%
%
%  Student's code starts here...
%
%
%
%  Good luck CE Student!!!
%
