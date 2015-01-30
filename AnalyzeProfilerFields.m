function varargout = AnalyzeProfilerFields(varargin)
% AnalyzeProfilerFields reads in a data structure (created via ParseSNCprm
% or ParseSNCtxt) and computes statistics for each field.  If multiple
% exposures are detected when analyzing time-dependent ParseSNCprm data
% (indicated by gaps in time between measurements), each profile will be
% separated out and statistics will be computed independently.  The
% various statistics are returned as a MATLAB structure.
%
% This function may be executed with a second input argument containing one
% or more reference profiles for each axis (xdata, ydata, pdiag, ndiag). If
% provided, additional statistics (such as Gamma comparisons) will be 
% computed. If only one reference profile is provided, but multiple
% measured profiles are provided, the same reference profile will be
% applied to each measured profile. If multiple reference profiles are
% provided, the closest matching reference profile to the measured profile 
% will be selected (using a correlation coefficient).
%
% If reference data is provided, a second structure can also be returned by
% this function containing the same fields (and array sizes) as the first
% results array, but containing the statistics for the reference profile.
%
% The following variables are required for proper execution:
%   varargin{1}: structure returned either by ParseSNCprm or ParseSNCtxt 
%       (see the documentation for these functions for more information on 
%       the fields contained)
%   varargin{2} (optional): structure containing reference profile data. 
%       The following structure fields can be incuded: xdata (2 x n+1 
%       array), ydata (2 x n+1 array), pdiag (2 x n+1 array), and ndiag 
%       (2 x n+1 array). To compute Gamma, must also include the fields
%       abs (absolute gamma criterion), dta (distance to agreement, in cm),
%       and local (flag indicating whether to perform local or global Gamma
%       comparison).
%
% The following structure fields are returned for each return variable upon 
% successful completion:
%   xdata: 2D array of background/ignore flag/calibration corrected X axis 
%       data, where column one is the position (in cm), and columns 2:n+1 
%       are the data for each field
%   ydata: 2D array of background/ignore flag/calibration corrected Y axis 
%       data, where column one is the position (in cm), and columns 2:n+1 
%       are the data for each field
%   pdiag: 2D array of background/ignore flag/calibration corrected 
%       positive diagonal axis data, where column one is the position (in 
%       cm), and columns 2:n+1 are the data for each field
%   ndiag: 2D array of background/ignore flag/calibration corrected 
%       negative diagonal axis data, where column one is the position (in 
%       cm), and columns 2:n+1 are the data for each field
%   tdata (optional): 2D array of time-dependent central channel data, 
%       where column one is the absolute time and column two is the 
%       differential central channel response. If time dependent data is
%       not available (i.e. ParseSNCtxt data), this field is not returned.
%   xfwhms: vector of X axis Full Width at Half Maximum(s) for each field
%   xedges: n x 2 array of left and right FWHM-defined X axis field edges
%   yfwhms: vector of Y axis Full Width at Half Maximum(s) for each field
%   yedges: n x 2 array of left and right FWHM-defined Y axis field edges
%   pfwhms: vector of positive diagonal axis Full Width at Half Maximum(s) 
%       for each field
%   pedges: n x 2 array of left and right FWHM-defined positive diagonal 
%       axis field edges
%   nfwhms: vector of negative diagonal axis Full Width at Half Maximum(s) 
%       for each field
%   nedges: n x 2 array of left and right FWHM-defined negative diagonal 
%       axis field edges
%   xflat: vector of central 80% flatness for X axis profiles
%   xsym: vector of areal symmetry for X axis profiles
%   yflat: vector of central 80% flatness for Y axis profiles
%   ysym: vector of areal symmetry for Y axis profiles
%   pflat: vector of central 80% flatness for positive diagonal profiles
%   psym: vector of areal symmetry for positive diagonal profiles
%   nflat: vector of central 80% flatness for negative diagonal profiles
%   nsym: vector of areal symmetry for negative diagonal profiles
%
% The following additional structure fields are returned in the first
% return variable if reference data was provided.
%   ref (optional): vector indicating which reference profile was selected
%   corr (optional): n x 4 array containing the correlation coefficient 
%       between the measured profiles and reference profile, for each axis
%       (x, y, p, and n)
%   xdiff (optional): 2D array of X axis differences, where column one 
%       is the position (in cm), and columns 2:n+1 are the abs differences
%   ydiff (optional): 2D array of Y axis differences, where column one 
%       is the position (in cm), and columns 2:n+1 are the abs differences
%   pdiff (optional): 2D array of positive diagonal differences, where 
%       column one is the position (in cm), and columns 2:n+1 are the abs
%       differences
%   ndiff (optional): 2D array of negative diagonal differences, where 
%       column one is the position (in cm), and columns 2:n+1 are the abs
%       differences
%   xgamma (optional): 2D array of X axis Gamma values, where column one 
%       is the position (in cm), and columns 2:n+1 are the Gamma indices
%   ygamma (optional): 2D array of Y axis Gamma values, where column one 
%       is the position (in cm), and columns 2:n+1 are the Gamma indices
%   pgamma (optional): 2D array of positive diagonal Gamma values, where 
%       column one is the position (in cm), and columns 2:n+1 are the Gamma
%       indices
%   ngamma (optional): 2D array of negative diagonal Gamma values, where 
%       column one is the position (in cm), and columns 2:n+1 are the Gamma
%       indices





