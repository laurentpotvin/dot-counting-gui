function [ offset_x, offset_y ] = align_frames( frame1, frame2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% raw_data = sum(frame1,3);
% raw_data(:,:,2)=sum(frame2,3);

raw_data = frame1(:,:,round(size(frame1,3)/2));
raw_data(:,:,2)= frame2(:,:,round(size(frame2,3)/2));

ROI = imSelectROI(raw_data(:,:,1));
 % Compute correlation function
 [Gtime] = sticsTFM(raw_data(ROI.Xrange,ROI.Yrange,:),size(raw_data(ROI.Xrange,ROI.Yrange,:),3));
%   [Gtime] = sticsTFM(raw_data(:,:,:),size(raw_data(:,:,:),3));

  % TO DO : better way to estimate binning and PSF width
  if size(raw_data,1) >1500
      beam_size = 75;
  else
      beam_size = 35;
  end
  
% TO DO: autocrop before fitting to speed up process
% Gaussian fitting
   [coeffGtime(1:size(Gtime,3),:)] = gaussfit(Gtime/(mean(mean(mean(Gtime)))),'time',1,'n',beam_size);
    
   offset_x = coeffGtime(1,5) - coeffGtime(2,5);
   offset_y = coeffGtime(1,6) - coeffGtime(2,6);
   
end

function [timecorr] = sticsTFM(imgser,upperTauLimit)

% July 10, 2003
% David Kolin
% Updated Nov 21, 2006 to allow for cross correlation
% Calculates the full time correlation function given 3D array of image series
% Usage:
% [timecorr] = stics(imgser,upperTauLimit) 
% OR
% [crossCorr] = stics(imgser1,imgser2,upperTauLimit)
% OR
% [crossCorr autoCorr1 autoCorr2] = stics(imgser1,imgser2,upperTauLimit)
% where timecorr is the spatio-temporal correlation function of imgser,
% calculated up to time lag upperTauLimit
% In the second case, timecorr is the spatio-temporal cross-correlation
% function of imgser1 and imgser2
% Modified April 13/07 to calculate xcorr 1 vs 2 and then 2 vs 1 and
% average them
% Updated July 15, 2009 to correct for edge effect

useWaitBar = 1;

if useWaitBar
    set(gcbf,'pointer','watch');
    h = waitbar(0,'Calculating time correlation functions...');
end

% tau is the lag
% pair is the nth pair of a lag time

% if length(varargin)==2
%     imgser = varargin{1};
%     upperTauLimit = min(varargin{2},size(imgser,3));
    
    %timecorr = zeros(size(imgser,1),size(imgser,2),upperTauLimit);  % preallocates lagcorr matrix for storing raw time corr functions
    %SeriesMean = squeeze(mean(mean(imgser)));
    fftStack = double(zeros(size(imgser)));
    
    for i = 1:size(imgser,3)
        fftStack(:,:,i) = fft2(double(imgser(:,:,i)));
    end
    
    for tau = 0:upperTauLimit-1
        lagcorr = zeros(size(imgser,1),size(imgser,2),(size(imgser,3)-tau));
        for pair=1:(size(imgser,3)-tau)
            lagcorr(:,:,pair) = ifft2mod(fftStack(:,:,pair).*conj(fftStack(:,:,(pair+tau))),'symmetric');
        end
        timecorr(:,:,(tau+1)) = fftshift(mean(lagcorr,3));
        
%         % Checks for significance of global maximum
%          if (~correlationSignificance(timecorr(:,:,tau+1)))
%              timecorr = timecorr(:,:,1:(end-1)); % cut off the "bad" lag
%              break
%          end
        
        if useWaitBar
            if ishandle(h)
                waitbar((tau+1)/(upperTauLimit),h)
            else
                break
            end
        end
    end


if useWaitBar
    close(h)
end
set(gcbf,'pointer','arrow');
end

function x = ifft2mod(varargin)
%IFFT2 Two-dimensional inverse discrete Fourier transform.
%   Same as original ifft2, but without error checking.
%   IFFT2(F) returns the two-dimensional inverse Fourier transform of matrix
%   F.  If F is a vector, the result will have the same orientation.
%
%   IFFT2(F,MROWS,NCOLS) pads matrix F with zeros to size MROWS-by-NCOLS
%   before transforming.
%
%   IFFT2(..., 'symmetric') causes IFFT2 to treat F as conjugate symmetric
%   in two dimensions so that the output is purely real.  This option is
%   useful when F is not exactly conjugate symmetric merely because of
%   round-off error.  See the reference page for the specific mathematical
%   definition of this symmetry.
%
%   IFFT2(..., 'nonsymmetric') causes IFFT2 to make no assumptions about the
%   symmetry of F.
%
%   Class support for input F:
%      float: double, single
%
%   See also FFT, FFT2, FFTN, FFTSHIFT, FFTW, IFFT, IFFTN.

%   Copyright 1984-2006 The MathWorks, Inc. 
%   $Revision: 5.16.4.4 $  $Date: 2006/10/02 16:32:15 $

%error(nargchk(1, 4, nargin, 'struct'))

f = varargin{1};
m_in = size(f, 1);
n_in = size(f, 2);
num_inputs = numel(varargin);
symmetry = 'nonsymmetric';
if ischar(varargin{end})
    symmetry = varargin{end};
    num_inputs = num_inputs - 1;
end

if num_inputs == 1
    m_out = m_in;
    n_out = n_in;

elseif num_inputs == 2
    error('MATLAB:ifft2:invalidSyntax', ...
          'If you specify MROWS, you also have to specify NCOLS.')
    
else
    m_out = double(varargin{2});
    n_out = double(varargin{3});
end

if ~isa(f, 'float')
    f = double(f);
end

if (m_out ~= m_in) || (n_out ~= n_in)
    out_size = size(f);
    out_size(1) = m_out;
    out_size(2) = n_out;
    f2 = zeros(out_size, class(f));
    mm = min(m_out, m_in);
    nn = min(n_out, n_in);
    f2(1:mm, 1:nn, :) = f(1:mm, 1:nn, :);
    f = f2;
end

if ndims(f)==2
    x = ifftn(f, symmetry);
else
    x = ifft(ifft(f, [], 2), [], 1, symmetry);
end   


end

function [a] = gaussfit(corr,type,pixelsize,whitenoise,radius)

% Usage: a = gaussfit(corr,type,pixelsize,whitenoise);

%set(gcbf,'pointer','watch');

[X,Y] = meshgrid(-((size(corr,2)-1)/2)*pixelsize:pixelsize:((size(corr,2)-1)/2)*pixelsize,-((size(corr,1)-1)/2)*pixelsize:pixelsize:(size(corr,1)-1)/2*pixelsize);
grid = [X Y];

[Y0, X0] = find(ismember(corr,max(max(corr))),size(corr,3));
X0 = mod(X0,size(corr,2));

% Find X0 and Y0 are where remainder from mod was zero -- these are set to
%the "max" (ie size) of the corr
X0(ismember(X0,0)) = size(corr,2);
X0(ismember(Y0,0)) = size(corr,1);

% Sets curve fit options, and sets lower bounds for amplitude and beam
% radius to zero
lb = [0 0 -1 min(min(grid)) min(min(grid))];
ub = [];

if nargin < 5
    weights = ones(size(corr));
else
    weights = ones(size(corr));
    for i=1:size(corr,3)
        weights(:,:,i) = circle(size(corr,2),size(corr,1),X0(i), Y0(i),radius);
    end
end
    
% If there's whitenoise, 2 highest values (NB this might be more than
% two points!) in corr func are set to zero, and given no weight in the fit

if strcmp(whitenoise,'y')&&strcmp(type,'2d')
    for j=1:1
            i = find(ismember(corr(:,:,:),max(max(corr(:,:,:)))));
            %ZerChan = i;
            corr(i) = 0;
            weights(i) = 0;
    end
end

y0 = min(min(corr));
y0 = squeeze(y0);
g0 = squeeze(max(max(corr))) - y0;

% wguess = zeros(size(corr,3),1);
% for i=1:size(corr,3)
% [Wy, Wx] = find(ismember(abs((corr(:,:,i)/g0(i) - exp(-1))),min(min(abs(corr(:,:,i)/g0(i) - exp(-1))))));
% Wx = mod(Wx,size(corr,2));
% wguess(i) = mean(( (Wx - X0(i)).^2  + (Wy - Y0(i)).^2   ).^(1/2))*pixelsize;
% end
wguess = 0.4*ones(size(g0));

% Converts from matrix index to LOCATION in pixelsize units
for i=1:size(corr,3)
    X0(i) = X(1,X0(i));
end
for i=1:size(corr,3)
    Y0(i) = Y(Y0(i),1);
end

options = optimset('Display','off');
displayWaitbar = 'n';

if strcmp(displayWaitbar,'y')
    h = waitbar(0,'Fitting correlation functions...');
end

a = zeros(size(corr,3),6);

% Fits each corr func separately
switch lower(type)
    case '2d'
        initguess = [g0 wguess y0 X0 Y0];
        for i=1:size(corr,3)
            if strcmp(displayWaitbar,'y'); waitbar(i/size(corr,3),h); end
            %a0 = initguess(i,:);
            a0xy(1:2) = initguess(i,1:2);
            a0xy(3) = a0xy(2);
            a0xy(4:6) = initguess(i,3:5);
            a(i,:) = lsqcurvefit(@gauss2dwxy,a0xy,grid,corr(:,:,i).*weights(:,:,i),lb,ub,options,weights(:,:,i));
        end
    case 'time'
        initguess = [g0 wguess wguess y0 X0 Y0];
        for i=1:size(corr,3)
            if strcmp(displayWaitbar,'y'); waitbar(i/size(corr,3),h); end
            if i==1
                a0 = initguess(i,:);
            else
                a0 = a(i-1,:);
            end
            % sneak weights into end of grid matrix
            grid(:,(size(X,2)*2+1):(size(X,2)*3)) = weights(:,:,i);
            funlist = {1, @(a,grid) exp(-((grid(:,1:size(grid,2)/3)-a(2)).^2+(grid(:,size(grid,2)/3+1:2*size(grid,2)/3)-a(3)).^2)/(a(1)^2)) .* grid(:,2*size(grid,2)/3+1:end)  };
            NLPstart = [a0(2) a0(5) a0(6)];
            warning('off','MATLAB:rankDeficientMatrix');
            [INLP,ILP] = pleas(funlist,NLPstart,grid,corr(:,:,i),options);
            warning('on','MATLAB:rankDeficientMatrix');
            a(i,1) = ILP(2);
            a(i,2) = INLP(1);
            a(i,3) = INLP(1);
            a(i,4) = ILP(1);
            a(i,5) = INLP(2);
            a(i,6) = INLP(3);
            % "old" way -- all parameters determined in a nonlinear fit
            %[a(i,:),res(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = lsqcurvefit(@gauss2d,a0,grid,corr(:,:,i).*weights(:,:,i),lb,ub,curvefitoptions,weights(:,:,i));
        end
    case 'timeasym'
        initguess = [g0 wguess wguess y0 X0 Y0];
        for i=1:size(corr,3)
            if strcmp(displayWaitbar,'y'); waitbar(i/size(corr,3),h); end
            if i==1
                a0 = initguess(i,:);
            else
                a0 = a(i-1,:);
            end
            funlist = {1, @(a,grid) exp(-(    ((grid(:,1:size(grid,2)/2)-a(3))/a(1)).^2  +  ((grid(:,size(grid,2)/2+1:end)-a(4))/a(2)).^2 ) )}; % instead of 0:  grid(:,1:size(grid,2)/2).*grid(:,size(grid,2)/2+1:end)
            NLPstart = [a0(2) a0(3) a0(5) a0(6)];
            warning('off','MATLAB:rankDeficientMatrix');
            [INLP,ILP] = pleas(funlist,NLPstart,grid,corr(:,:,i),options);
            warning('on','MATLAB:rankDeficientMatrix');
            a(i,1) = ILP(2);
            a(i,2) = INLP(1);
            a(i,3) = INLP(2);
            a(i,4) = ILP(1);
            a(i,5) = INLP(3);
            a(i,6) = INLP(4);
            %[a(i,:),res(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = lsqcurvefit(@gauss2d,a0,grid,corr(:,:,i).*weights(:,:,i),lb,ub,curvefitoptions,weights(:,:,i));
        end
    otherwise
        error('Fitting mode must be ''2d'', ''time'', or ''timeasym''.');
end

% If the peak moves "past" the edge of the correlation function, it will
% appear on the other side; this unwraps the positions so there are not
% discontinuities in the Gaussian position.  Does it separately for the x-
% and y-coordinates.
if strcmpi(type,'time') || strcmpi(type,'timeasym')
    a(:,5) = unwrapCustom(a(:,5),size(corr,2)/2*pixelsize,1);
    a(:,6) = unwrapCustom(a(:,6),size(corr,1)/2*pixelsize,1);
end

try close(h); catch end % if it waitbar was closed or never open to begin with

%set(gcbf,'pointer','arrow');
end

function [c_mask]=circle(ix,iy,cx,cy,r)

% By M. Bach
% Draws a circle with centre cx,cy in image ix,iy with radius r
% [c_mask]=circle(ix,iy,cx,cy,r)

[x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
c_mask=((x.^2+y.^2)<=r^2);

end
function [INLP,ILP] = pleas(funlist,NLPstart,xdata,ydata,options)
% pleas: partitioned nonlinear least squares estimation
% usage: [INLP,ILP] = pleas(funlist,NLPstart,xdata,ydata)
% usage: [INLP,ILP] = pleas(funlist,NLPstart,xdata,ydata,options)
%
% arguments: (input)
%  funlist - cell array of functions comprising the nonlinear parts
%            of each term in the model. Each independent function in
%            this list must transform xdata using a vector of intrinsicly
%            nonlinear parameters into an array of the same size and
%            shape as ydata. The arguments to each function will be in
%            the order (coef,xdata).
%
%            These functions may be
%             - scalar (double) constants (E.g., 1)
%             - anonymous functions
%             - inline functions
%             - character function names
%
%  NLPstart - vector of starting values for the intrinsicly nonlinear
%            parameters only.
%
%  xdata   - array of independent variables
%            
%  ydata   - array of dependent variable data
%
%  options - options structure appropriate for lsqnonlin
%
%
% arguments (output)
%  INLP - optimized list of intrinsicly nonlinear parameters
%
%  ILP  - optimized list of intrinsicly linear parameters 
%
%
% Example usage:
%  Fit a simple exponential model plus a constant term to data
%
%   x = rand(100,1);
%   y = 4 - 3*exp(2*x) + randn(size(x));
%
%   funlist = {1, @(xdata,coef) exp(xdata*coef)};
%   NLPstart = 1;
%   options = optimset('disp','iter');
%   [INLP,ILP] = pleas(funlist,NLPstart,x,y,options)
%
% Output:
%                                          Norm of      First-order 
%  Iteration  Func-count     f(x)          step          optimality   CG-iterations
%     0          2         116.796                          40.5
%     1          4         74.6378        1.00406            2.6            1
%     2          6         74.4382      0.0758513         0.0443            1
%     3          8         74.4381     0.00131311       0.000597            1
% Optimization terminated: relative function value
%  changing by less than OPTIONS.TolFun.
%
% INLP =
%    2.0812
%
% ILP =
%    3.6687
%   -2.7327

% Check the functions
if ~iscell(funlist)
  error 'funlist must be a cell array of functions, even if only one fun'
end
nfun=length(funlist);
for i = 1:nfun
  fi = funlist{i};
  
  % There are two cases where we need to turn the supplied
  % function into an executable function
  if isa(fi,'double')
    % a constant
    funlist{i} = @(xdata,coef) repmat(fi,size(ydata));
  elseif ischar(fi)
    % a character function name
    funlist{i} = str2func(fi);
  end
end

% were any options supplied?
if (nargin<5) || isempty(options)
  options = optimset('lsqnonlin');
end

% make sure that ydata is a column vector
ydata = ydata(:);
ny = length(ydata);

% ================================================
% =========== begin nested function ==============
% ================================================
function [res,ILP] = pleas_obj(INLP)
  % nested objective function for lsqnonlin, so all
  % the data and funs from pleas are visible to pleas_obj
  
  % loop over funlist
  A = zeros(ny,nfun);
  for i=1:nfun
    fi = funlist{i};
    term = fi(INLP,xdata);
    A(:,i) = term(:);
  end
  
  % do the linear regression using \
  ILP = A\ydata;
  
  % residuals for lsqnonlin
  res = A*ILP - ydata;
  
end % nested function termination
% ================================================
% ============= end nested function ==============
% ================================================

% call lsqnonlin, using a nested function handle
LB=[];
UB=[];
INLP = lsqnonlin(@pleas_obj,NLPstart,LB,UB,options);

% call one final time to get the final linear parameters
[junk,ILP] = pleas_obj(INLP);

end % main function terminator

function q = unwrapCustom(p,cutoff,dim)
%UNWRAPCUSTOM Unwrap phase angle.
%   UNWRAPCUSTOM(P) unwraps radian phases P by changing absolute 
%   jumps greater than or equal to pi to their 2*pi complement.
%   It unwraps along the first non-singleton dimension of P.
%   P can be a scalar, vector, matrix, or N-D array. 
%
%   UNWRAPCUSTOM(P,TOL) uses a jump tolerance of TOL rather
%   than the default TOL = pi.
%
%   UNWRAPCUSTOM(P,[],DIM) unwraps along dimension DIM using the
%   default tolerance. UNWRAP(P,TOL,DIM) uses a jump tolerance
%   of TOL.
%
%   Unlike UNWRAP, UNWRAPCUSTOM will unwrap for values smaller than pi, and
%   will change jumps larger than P to their 2*P complement.
%
%   Class support for input P:
%      float: double, single
%
%   See also ANGLE, ABS.

%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 5.14.4.3 $  $Date: 2005/11/18 14:14:38 $
%   Modified July 17, 2007 by David Kolin

% Overview of the algorithm:
%    Reshape p to be a matrix of column vectors. Perform the 
%    unwrap calculation column-wise on this matrix. (Note that this is
%    equivalent to performing the calculation on dimension one.) 
%    Then reshape the output back.

ni = nargin;

% Treat row vector as a column vector (unless DIM is specified)
rflag = 0;
if ni<3 && (ndims(p)==2) && (size(p,1)==1), 
   rflag = 1; 
   p = p.';
end

% Initialize parameters.
nshifts = 0;
perm = 1:ndims(p);
switch ni
case 1
   [p,nshifts] = shiftdim(p);
   cutoff = pi;     % Original UNWRAP used pi*170/180.
case 2
   [p,nshifts] = shiftdim(p);
otherwise    % nargin == 3
   perm = [dim:max(ndims(p),dim) 1:dim-1];
   p = permute(p,perm);
   if isempty(cutoff),
      cutoff = pi; 
   end
end
   
% Reshape p to a matrix.
siz = size(p);
p = reshape(p, [siz(1) prod(siz(2:end))]);

% Unwrap each column of p
q = p;
for j=1:size(p,2)
   % Find NaN's and Inf's
   indf = find(isfinite(p(:,j)));
   % Unwrap finite data (skip non finite entries)
   q(indf,j) = LocalUnwrap(p(indf,j),cutoff);
end

% Reshape output
q = reshape(q,siz);
q = ipermute(q,perm);
q = shiftdim(q,-nshifts);
if rflag, 
   q = q.'; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local Functions  %%%%%%%%%%%%%%%%%%%

function p = LocalUnwrap(p,cutoff)
%LocalUnwrap   Unwraps column vector of phase values.

m = length(p);

% Unwrap phase angles.  Algorithm minimizes the incremental phase variation 
% by constraining it to the range [-pi,pi]
dp = diff(p,1,1);                        % Incremental phase variations
dps = mod(dp+cutoff,2*cutoff) - cutoff;  % Equivalent phase variations in [-pi,pi)
dps(dps==-cutoff & dp>0,:) = cutoff;     % Preserve variation sign for pi vs. -pi
dp_corr = dps - dp;                      % Incremental phase corrections
dp_corr(abs(dp)<cutoff,:) = 0;           % Ignore correction when incr. variation is < CUTOFF

% Integrate corrections and add to P to produce smoothed phase values
p(2:m,:) = p(2:m,:) + cumsum(dp_corr,1);
end
end


function varargout = imSelectROI(img,varargin)
%==========================================================================
%  Author: Andriy Nych ( nych.andriy@gmail.com )
% Version: 733250.03649467591
%--------------------------------------------------------------------------
% This functions displays GUI for selecting square or rectangular part
% of the input image IMG. To perform selection user must click mouse twice:
% at two corners of the selection area.
% User can change the shape at any moment, even when first point is set,
% unless it is not forbidden by additional parameters.
% Use also cam change the way the selection area is calculated
% from the two selected points.
% Depending on the combination of the shape and mode it could be:
%--------------------------------------------------------------------------
% Shape       Mode        Result
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Rectangle   Free        Rectangle with one corner at first point (P1)
%                         and another corner at second point (P2).
% Rectangle   Centered    Rectangle with its center at first point (P1)
%                         and corner at second point (P2).
% Square      Free        Square of the largest size that can be
%                         fitted into rectangle made by (P1) and (P2)
%                         with one corner at (P1).
% Square      Centered    Square of the largest size that can be
%                         fitted into centered rectangle.
%                         Center of the square is at (P1).
%--------------------------------------------------------------------------
% Behavior of the imSelectROI can be altered by providing additional
% parameters in MatLab well-known ParamName,ParamValue style.
%
% NOTE      This function was developed under MatLab R2006b.
% ====      It requires Image Processing Toolbox to work.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Syntax
% ======
%                   imSelectROI( img, param1, val1, param2, val2, ...)
%             ROI = imSelectROI( img, param1, val1, param2, val2, ...)
% [ ROI, SelImg ] = imSelectROI( img, param1, val1, param2, val2, ...)
%
% Displays GUI and returns:
%
% SelImg - Selected part of the image passed as first parameter.
%          Even if first parameter is a file name (see below).
%
% ROI - structure with fields:
%   ROI.Xmin    - minimal value of X coordinate of the selected area
%   ROI.Xmax    - maximal value of X coordinate of the selected area
%   ROI.Ymin    - minimal value of Y coordinate of the selected area
%   ROI.Ymax    - maximal value of Y coordinate of the selected area
%   ROI.DX      - horizontal size of the selected area
%       ROI.DX = ROI.Xmax - ROI.Xmin + 1
%   ROI.DY      - vertical size of the selected area
%       ROI.DY = ROI.Ymax - ROI.Ymin + 1
%   ROI.Xrange  - same as [ROI.Xmin:ROI.Xmax]
%   ROI.Yrange  - same as [ROI.Ymin:ROI.Ymax]
%
%   Selected part can be retrieved from original image as
%       img( ROI.Xrange, ROI.Yrange, :)
%   This allows to perform selection once and use the same ROI
%   to process series of images (see examples at hte end).
%
% Arguments
% =========
%
% img     Anything that can be passed to IMSHOW as a single parameter.
%         In could be file name or preloaded image.
%         See "help imshow" for more information about the syntaxes.
%
% Parameters
% ==========
%
% AllowedShape  (string): {'Any'} | 'Square' | 'Rectangle'
%
%   This parameter controls shape of the selection.
%   Specifying 'Square' or 'Rectangle' you prevent user from
%   selecting other shape.
%   By specifying 'Any' or omitting 'AllowedShape' at all
%   user is allowed to select any shape.
%
% SelectionMode (string): {'Free'} | 'Centered'
%
%   This parameter controls selection mode.
%   But in this case user still can select other mode.
%
% FastReturn    (string): {'off'} | 'on'
%
%   This parameter controls how the GUI behaves when user finishes
%   seletion. 
%   When 'off' value provided function waits for user to press
%   "DONE" button, allowing user to change selection by
%   "START OVER" button.
%   When 'on' value provided function returns immediately after user
%   makes valid selection of second point. In this case it is also
%   possible to change selection, but only until the second point was
%   NOT selected by user.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Examples
% ======
% ROI = imSelectROI( 'c:\Image.jpg');
% [ROI,SelImage] = imSelectROI( 'c:\Image.jpg', 'AllowedShape','Square');
%
% % FNames is a cell array of image file names
% ROI = imSelectROI( FNames{1} );
% for i=1:length(FNames)
%     image = imread(FNames{i}); %whole image
%     selection = image( ROI.Xrange, ROI.Yrange, :); %selected area
%     ...
% end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Additional info
% ===============
%   help imshow     for additional information on what can be passed
%                   to imSelectROI as first argument.
%==========================================================================

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% GUI parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LBoxWidth       = 0.1;  % width of the ListBoxes in NORMALIZED units
                        % Change this when width of listboxes is too smal
                        % to accomodate its strings or 
                        % to meet your preferences
LBoxHeight      = 44;   % height of the ListBoxes in PIXELS
ScrMargin       = 0.05; % Screen margin for the whole figure
                        % in NORMALIZED units

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parsing arguments
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% GP.AllowedShapes    = {'Any','Square','Rectangle'};
GP.AllowedShapes    = {'Square','Rectangle'};
GP.SelectionModes   = {'Free','Centered'};
GP.FastReturns      = {'on','off'};

GP.AllowedShape_i   = GetParameterValue('AllowedShape',varargin, 'Any');
GP.SelectionMode_i  = GetParameterValue('SelectionMode',varargin, 'Free');
GP.FastReturn_i     = GetParameterValue('FastReturn',varargin, 'Off');

if any( strcmpi( GP.AllowedShape_i,  {GP.AllowedShapes{:} 'Any'} ) )
    GP.SelectionMode    = GP.SelectionMode_i;
else
    error('%s: Wrong value of "AllowedShape" parameter ("%s")',mfilename,GP.AllowedShape_i);
end

if any( strcmpi( GP.SelectionMode_i, GP.SelectionModes ) )
    GP.AllowedShape     = GP.AllowedShape_i;
else
    error('%s: Wrong value of "SelectionMode" parameter ("%s")',mfilename,GP.SelectionMode_i);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Creating GUI
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GP.hFig         = figure(...
    'Name','Select subimage',...
    'NumberTitle','off',...
    'Color',[1 1 1],...
    'WindowStyle','modal',...{normal} | 
    'Pointer','crosshair',... fullcrosshair crosshair
    'Resize','off',...{on} | off
    'Toolbar','none',...
    'Menubar','none');
GP.hImg         = imagesc(img);
GP.XLim         = xlim;
GP.YLim         = ylim;
GP.hAxes        = gca;

set(GP.hAxes, 'ActivePositionProperty','OuterPosition');

GP.hPixVal      = impixelinfoval(GP.hFig,GP.hImg);

set(GP.hFig,...
    'Units','normalized',...
    'Position',[ ScrMargin ScrMargin 1-2*ScrMargin 1-2*ScrMargin ]);

GP.hSP          = imscrollpanel(GP.hFig,GP.hImg);

set(GP.hSP,...
    'Units','normalized',...
    'Position',[0 .1 1 .9]);
GP.hSPAPI       = iptgetapi(GP.hSP);
GP.hSPAPI.setMagnification(1);

movegui(GP.hFig, 'east');

GP.hDoneBtn     = uicontrol( 'Style','pushbutton',  'Parent',gcf,   'FontWeight','bold',                'units','normalized',   'position',[1-2.0*LBoxWidth 0 LBoxWidth 0.2], 'string','Done');
GP.hModeLBox    = uicontrol( 'Style','listbox',     'Parent',gcf,   'min',0,    'max',1,    'value',1,  'units','normalized',   'position',[1-3.5*LBoxWidth 0 LBoxWidth 0.2], 'string',GP.SelectionModes);
GP.hShapeLBox   = uicontrol( 'Style','listbox',     'Parent',gcf,   'min',0,    'max',1,    'value',1,  'units','normalized',   'position',[1-4.5*LBoxWidth 0 LBoxWidth 0.2], 'string',GP.AllowedShapes);
GP.hZoomPMenu   = uicontrol( 'Style','popupmenu',   'Parent',gcf,                           'value',5,  'units','normalized',   'position',[1-6.0*LBoxWidth 0 LBoxWidth 0.2], 'string',{'0010%','0025%','0050%','0075%','0100%','0150%','0200%','0250%','0300%','0400%','0500%','0600%','0700%','0800%','0900%','1000%'});
GP.hRedoBtn     = uicontrol( 'Style','pushbutton',  'Parent',gcf,   'FontWeight','bold',                'units','normalized',   'position',[1-7.5*LBoxWidth 0 LBoxWidth 0.2], 'string','Start over');
GP.hStatus      = uicontrol( 'Style','text',        'Parent',gcf,   'FontWeight','bold', 'FontSize',12, 'units','normalized',   'position',[1-9.0*LBoxWidth 0 LBoxWidth 0.2], 'string','');

set(GP.hShapeLBox,  'units','pixels'); p = get(GP.hShapeLBox, 'position'); set(GP.hShapeLBox, 'position',[ p(1)  0               p(3)    LBoxHeight ]);
set(GP.hModeLBox,   'units','pixels'); p = get(GP.hModeLBox,  'position'); set(GP.hModeLBox,  'position',[ p(1)  0               p(3)    LBoxHeight ]);
set(GP.hDoneBtn,    'units','pixels'); p = get(GP.hDoneBtn,   'position'); set(GP.hDoneBtn,   'position',[ p(1)  0               p(3)    LBoxHeight ]);
set(GP.hZoomPMenu,  'units','pixels'); p = get(GP.hZoomPMenu, 'position'); set(GP.hZoomPMenu, 'position',[ p(1)  LBoxHeight-25   p(3)    25         ]);
set(GP.hRedoBtn,    'units','pixels'); p = get(GP.hRedoBtn,   'position'); set(GP.hRedoBtn,   'position',[ p(1)  0               p(3)    LBoxHeight ]);
set(GP.hStatus,     'units','pixels'); q = get(GP.hStatus,    'position'); set(GP.hStatus,    'position',[ q(1)  LBoxHeight+08   p(3)*8  25         ]);

if strcmpi(GP.AllowedShape_i,'Square')
    % Square shape
    set(GP.hShapeLBox, 'Value',1);
    set(GP.hShapeLBox, 'enable','inactive');
    set(GP.hShapeLBox, 'TooltipString','Only square selection allowed');
elseif strcmpi(GP.AllowedShape_i,'Rectangle')
    % Rectangle shape
    set(GP.hShapeLBox, 'Value',2);
    set(GP.hShapeLBox, 'enable','inactive');
    set(GP.hShapeLBox, 'TooltipString','Only rectangular selection allowed');
else
    % Any shape
    set(GP.hShapeLBox, 'Value',2);
    set(GP.hShapeLBox, 'TooltipString','Select shape of your selection');
end

if strcmpi(GP.SelectionMode_i,'Centered')
    set(GP.hModeLBox, 'Value',2);
end

axes(GP.hAxes);
axis image; axis on;

%=========================================================================
% Stroing handles and other data to safe place
%=========================================================================
GP.P1   = [];
GP.P2   = [];
GP.SP1  = [];
GP.SP2  = [];
GP.SW   = [];
GP.SH   = [];
GP.hCrs = [];
GP.hRct = [];
guidata(gcf,GP);

%=========================================================================
% Tuning interface
%=========================================================================
set( gcf,           'WindowButtonUpFcn',    @FigureMouseBtnUp);
set( GP.hZoomPMenu, 'CallBack',             @AdjustZoom);
set( GP.hDoneBtn,   'CallBack',             @SelectionDone);
set( GP.hRedoBtn,   'CallBack',             @StartOver);

set( GP.hShapeLBox, 'CallBack',             @UpdateShapeAndMode);
set( GP.hModeLBox,  'CallBack',             @UpdateShapeAndMode);

iptaddcallback(gcf, 'WindowButtonMotionFcn', @FigureMouseMove);

uiwait;

if ishandle(GP.hFig)
    GP = guidata(gcf);
    delete(gcf);
%     ROI.Xmin    = GP.SP1(1);
%     ROI.Xmax    = GP.SP2(1);
%     ROI.Ymin    = GP.SP1(2);
%     ROI.Ymax    = GP.SP2(2);
    ROI.Xmin    = GP.SP1(2);
    ROI.Xmax    = GP.SP2(2);
    ROI.Ymin    = GP.SP1(1);
    ROI.Ymax    = GP.SP2(1);

    ROI.DX      = GP.SW;
    ROI.DY      = GP.SH;
    % fprintf(' nargin = %g\n',nargin);
    % fprintf('nargout = %g\n',nargout);
    switch nargout
        case 0, % No output parameters
            disp(ROI);
            if ischar(img)
                % FileName was provided
                im = imread(img);
                figure;
                imshow( im( ROI.Xmin:ROI.Xmax, ROI.Ymin:ROI.Ymax, : ) );
            else
                % Image matrix was provided
                figure;
                imshow( img( ROI.Xmin:ROI.Xmax, ROI.Ymin:ROI.Ymax, : ) );
            end
            return
        case 1, % 1 output parameter
            ROI.Xrange  = ROI.Xmin:ROI.Xmax;
            ROI.Yrange  = ROI.Ymin:ROI.Ymax;
            varargout{1} = ROI;
        case 2, % 2 output parameters
            ROI.Xrange  = ROI.Xmin:ROI.Xmax;
            ROI.Yrange  = ROI.Ymin:ROI.Ymax;
            varargout{1} = ROI;
            if ischar(img)
                % FileName was provided
                im = imread(img);
                varargout{2} = im( ROI.Xmin:ROI.Xmax, ROI.Ymin:ROI.Ymax, : );
            else
                % Image matrix was provided
                varargout{2} = img( ROI.Xmin:ROI.Xmax, ROI.Ymin:ROI.Ymax, : );
            end
        otherwise
            msgId = 'VarArgOutTest:invalidNumOutputArguments';
            msg = 'Tne number of output arguments is invalid.';
            error(msgId,'%s',msg);
    end

else
    fprintf('%s: Warning! Figure was closed before selection was retrieved! Empty matrix returned.\n',mfilename);
    switch nargout
        case 0,
            return;
        case 1,
            varargout{1} = [];
            return;
        case 2,
            varargout{1} = [];
            varargout{2} = [];
            return;
    end
end

return;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Done button callback
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function SelectionDone(src,evt)
PDelay  = 0.1;
GP = guidata(gcf);
if ~ValidateSelection
    tf = isempty( get(GP.hStatus, 'String') );
    if tf,set(GP.hStatus, 'String','ERROR: Wrong selection!'); end
    for i=1:7
        set(GP.hStatus, 'ForegroundColor',[1 0 0]); pause(PDelay);
        set(GP.hStatus, 'ForegroundColor',[0 0 0]); pause(PDelay);
    end
    if tf, set(GP.hStatus, 'String',''); end
else
    uiresume;
end
return;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% StartOver button callback
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function StartOver(src,evt)
GP = guidata(gcf);
GP.P1 = [];
GP.P2 = [];
GP.P1   = [];
GP.P2   = [];
GP.SP1  = [];
GP.SP2  = [];
GP.SW   = [];
GP.SH   = [];
set(gcf, 'Name',sprintf('[ %06g : %06g ]; P1 = [ %s]; P2 = [ %s]; Width = %06g; Height = %06g;',GP.cx,GP.cy, sprintf('%06g ',GP.P1), sprintf('%06g ',GP.P2), GP.SW, GP.SH)  );
set(GP.hStatus, 'String','');
DeleteHandleSafely(GP.hCrs); GP.hCrs = [];
DeleteHandleSafely(GP.hRct); GP.hRct = [];
guidata(gcf,GP);
return;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Mouse click callback
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function FigureMouseBtnUp(src,evt)
GP = guidata(gcf);
if isnan(GP.cx)||isnan(GP.cy)
    set(gcf, 'Name','MouseUp outside axes' );
else
    %set(gcf, 'Name',sprintf('[ %06g : %06g ]; P1 = [ %s]; P2 = [ %s]; Width = %06g; Height = %06g;',GP.cx,GP.cy, sprintf('%06g ',GP.P1), sprintf('%06g ',GP.P2), GP.SW, GP.SH)  );
    if isempty(GP.P1)
        % P1 is empty
        GP.P1 = [ GP.cx GP.cy ];
        set(gcf, 'Name',sprintf('[ %06g : %06g ]; P1 = [ %s]; P2 = [ %s]; Width = %06g; Height = %06g;',GP.cx,GP.cy, sprintf('%06g ',GP.P1), sprintf('%06g ',GP.P2), GP.SW, GP.SH)  );
        DeleteHandleSafely(GP.hCrs); GP.hCrs = [];
        if strcmpi( GP.SelectionMode, 'Centered' )
            GP.hCrs = caDrawCross(...
                GP.P1(1),GP.P1(2),13,...
                'xor', 1, 'k',':',...
                'k',8,'normal',...
                '','','','');
        end
        guidata(gcf,GP);
        return;
    else
        % P1 is NOT empty
        if isempty(GP.P2)
            % P2 is empty
            GP.P2 = [ GP.cx GP.cy ];
            %set(gcf, 'Name',sprintf('[ %06g : %06g ] P1 = [ %s] P2 = [ %s]',GP.cx,GP.cy, sprintf('%06g ',GP.P1), sprintf('%06g ',GP.P2) )  );
            %set(gcf, 'Name',sprintf('[ %06g : %06g ]; P1 = [ %s]; P2 = [ %s]; Width = %06g; Height = %06g;',GP.cx,GP.cy, sprintf('%06g ',GP.P1), sprintf('%06g ',GP.P2), GP.SW, GP.SH)  );
            %set(gcf, 'Name',sprintf('[ %06g : %06g ]; P1 = [ %s]; P2 = [ %s]; Width = %06g; Height = %06g;',GP.cx,GP.cy, sprintf('%06g ',GP.SP1), sprintf('%06g ',GP.SP2), GP.SW, GP.SH)  );
            set(gcf, 'Name',sprintf('P1 = [ %s]; P2 = [ %s]; Width = %06g; Height = %06g;',sprintf('%06g ',GP.SP1), sprintf('%06g ',GP.SP2), GP.SW, GP.SH)  );
            DeleteHandleSafely(GP.hCrs); GP.hCrs = [];
            DeleteHandleSafely(GP.hRct); GP.hRct = [];
            GP.hRct = caDrawRectangle(...
                GP.SP1(1),GP.SP1(2),...
                GP.SP2(1),GP.SP2(2),...
                'xor', 1, 'k', '-', 'none',...
                'k',8,'normal',...
                'selection','','','');
            guidata(gcf,GP);
            % FastReturn trick
            if strcmpi(GP.FastReturn_i,'on')
                % GP.FastReturn_i
                SelectionDone(src,evt);
            end
        else
            % P2 is NOT empty
            % set(gcf, 'Name',sprintf('P1 = [ %s]; P2 = [ %s]; Width = %06g; Height = %06g;',sprintf('%06g ',GP.SP1), sprintf('%06g ',GP.SP2), GP.SW, GP.SH)  );
        end
    end
end
return;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Mouse move callback
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function FigureMouseMove(src,evt)
GP = guidata(gcf);
% set(gcf, 'Name',get(GP.hPixVal,'String'));
ts = get(GP.hPixVal,'String');
p1 = strfind(ts, '(');
p2 = strfind(ts, ',');
p3 = strfind(ts, ')');
GP.cx   = str2double( ts( p1(1)+1 : p2(1)-1 ) );
GP.cy   = str2double( ts( p2(1)+1 : p3(1)-1 ) );
GP.CP   = [ GP.cx GP.cy ];

if ~isempty( GP.P1 )
    if ~isempty( GP.P2 )
        % P1 - P2
        DeleteHandleSafely(GP.hCrs); GP.hCrs = [];
        %set(gcf, 'Name',sprintf('P1 = [ %s]; P2 = [ %s]; Width = %06g; Height = %06g;', sprintf('%06g ',GP.SP1), sprintf('%06g ',GP.SP2), GP.SW, GP.SH)  );
    else
        % P1 - CP
        if strcmpi( GP.SelectionMode, 'Centered' )
            if isempty(GP.hCrs)
                GP.hCrs = caDrawCross(...
                    GP.P1(1),GP.P1(2),13,...
                    'xor', 1, 'k',':',...
                    'k',8,'normal',...
                    '','','','');
            end
        else
            DeleteHandleSafely(GP.hCrs); GP.hCrs = [];
        end
        set(GP.hStatus, 'String','');
        if isnan(GP.cx)||isnan(GP.cy)
            set(GP.hStatus, 'String','WARNING: Selection is partially outside the image');
        else
            [GP.SP1,GP.SP2] = RecalculateSelection(GP.P1,GP.CP, GP.AllowedShape, GP.SelectionMode);
            GP.SW = GP.SP2(1) - GP.SP1(1) + 1;
            GP.SH = GP.SP2(2) - GP.SP1(2) + 1;
            guidata(gcf,GP);
            if ~ValidateSelection
                set(GP.hStatus, 'String','WARNING: Selection is partially outside the image');
            end
            DeleteHandleSafely(GP.hRct); GP.hRct = [];
            %caDrawRectangle(X1,Y1,X2,Y2, EraseMode, LWidth,LColor,LStyle, FaceColor, TColor,TSize,FontWeight, TopLabel,BottomLabel,RightLabel,LeftLabel)
            GP.hRct = caDrawRectangle(...
                GP.SP1(1),GP.SP1(2),...
                GP.SP2(1),GP.SP2(2),...
                'xor', 1, 'k', ':', 'none',...
                'k',8,'normal',...
                'selection','','','');
            guidata(gcf,GP);
        end
        set(gcf, 'Name',sprintf('[ %06g : %06g ]; P1 = [ %s]; P2 = [ %s]; Width = %06g; Height = %06g;',GP.cx,GP.cy, sprintf('%06g ',GP.P1), sprintf('%06g ',GP.P2), GP.SW, GP.SH)  );
    end
    %set(gcf, 'Name',sprintf('[ %06g : %06g ]; P1 = [ %s]; P2 = [ %s]; Width = %06g; Height = %06g;',GP.cx,GP.cy, sprintf('%06g ',GP.P1), sprintf('%06g ',GP.P2), GP.SW, GP.SH)  );
end
guidata(gcf,GP);
return;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Zoom callback
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function AdjustZoom(src,evt)
GP  = guidata(gcf);
c   = get(gcbo,'String');
GP.hSPAPI.setMagnification( str2double( c{get(gcbo,'Value')}(1:end-1) ) / 100 );
return;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Deleting handle safely
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function DeleteHandleSafely(h)
for i=1:length(h(:))
    if ishandle( h(i) )
        delete ( h(i) );
    end
end
return;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Updating AllowedShape and SelectionMode
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function UpdateShapeAndMode(src,evt)
GP = guidata(gcf);
GP.AllowedShape     = GP.AllowedShapes  { get(GP.hShapeLBox,'Value') };
GP.SelectionMode    = GP.SelectionModes { get(GP.hModeLBox, 'Value') };
guidata(gcf,GP);
return;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Retrieving parameters from varargin
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function res = GetParameterValue(PName,PNameValArray,DefaultValue)
res = DefaultValue;   % in case PName is not present in PNameValArray
for i=1:2:length(PNameValArray)
    if strcmpi(PNameValArray{i},PName)
        if length(PNameValArray)>i
            res = PNameValArray{i+1};
            break;
        else
            error('%s: Parameter "%s" present, but its value absent',mfilename, PName);
        end
    end
end
return;
end
%==========================================================================
% Validate selection
%==========================================================================
function res = ValidateSelection
GP = guidata(gcf);
res = true;
if ~isfield(GP,'SP1')
    res = false; return;
end
if ~isfield(GP,'SP2')
    res = false; return;
end
res = res & ( GP.SP1(1) >= GP.XLim(1) );
res = res & ( GP.SP1(1) <= GP.XLim(2) );
res = res & ( GP.SP1(2) >= GP.YLim(1) );
res = res & ( GP.SP1(2) <= GP.YLim(2) );
res = res & ( GP.SP2(1) >= GP.XLim(1) );
res = res & ( GP.SP2(1) <= GP.XLim(2) );
res = res & ( GP.SP2(2) >= GP.YLim(1) );
res = res & ( GP.SP2(2) <= GP.YLim(2) );
return;
end
%==========================================================================
% Recalculate control points
%==========================================================================
function [r1,r2] = RecalculateControlPoints(p1,p2)
% fprintf( 'p1=[ %06g %06g ] p2=[ %06g %06g ]\n',p1,p2 );
r1  = [ min([ p1(1) p2(1) ]) min([ p1(2) p2(2) ]) ];
r2  = [ max([ p1(1) p2(1) ]) max([ p1(2) p2(2) ]) ];
% fprintf( 'r1=[ %06g %06g ] r2=[ %06g %06g ]\n',r1,r2 );
return;
end
%==========================================================================
% Recalculate selection rectangle
%==========================================================================
function [sp1,sp2] = RecalculateSelection(p1,p2, SelShape, SelMode)
XLen    = abs( p2(1) - p1(1) );
YLen    = abs( p2(2) - p1(2) );
MinLen  = min( [ XLen YLen ] );
Xsgn    = sign( p2(1) - p1(1) );
Ysgn    = sign( p2(2) - p1(2) );
if strcmpi(SelMode,'Centered')
    if strcmpi(SelShape,'Square')
        sp  = p1 - [Xsgn Ysgn]*MinLen;
        ep  = p1 + [Xsgn Ysgn]*MinLen;
    else
        sp  = p1 + ( p1 - p2 );
        ep  = p2;
    end
else
    if strcmpi(SelShape,'Square')
        sp  = p1;
        ep  = p1 + [Xsgn Ysgn]*MinLen;
    else
        sp  = p1;
        ep  = p2;
    end
end
[sp1,sp2] = RecalculateControlPoints(sp,ep);
return;
end
%==========================================================================
% Drawing cross in current axes
%==========================================================================
function h = caDrawCross(x,y,L, EraseMode, LWidth,LColor,LStyle, TColor,TSize,FontWeight, TLLabel,TRLabel,BLLabel,BRLabel)
%==========================================================================
%   This functions draws cross with four labels in current axes
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% h = caDrawCross(x,y,L, EraseMode, LWidth,LColor,LStyle,...
%                       TColor,TSize,FontWeight, TLLabel,TRLabel,BLLabel,BRLabel)
% 
%       x,y             Coordinates of center of cross
%       L               Cross size
%                           If L=0 this function draws cross trough full
%                           range of the coordinate system (like ginput)
%       EraseMode       EraseMode property of lines ad all text labels
%       LWidth          Line width
%       LColor          Line color
%       LStyle          Line style
%       TColor          Text color
%       TSize           Text size
%       FontWeight      Font weight
%       TLLabel         Top-left label string
%       TRLabel         Top-right label string
%       BLLabel         Bottom-left label string
%       BRLabel         Bottom-right label string
% 
%       h               array of handles of all created elements
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Type
%           help line
%   or
%           help text
%   for more information about line and text parameters
%==========================================================================
c   = 1;
if L>0
    h(c)    = line( [x-0 x+0],[y-L y+L+1], 'LineWidth',LWidth, 'LineStyle',LStyle, 'Color',LColor, 'EraseMode',EraseMode );
    c = c + 1;
    h(c)    = line( [x-L x+L+1],[y-0 y+0], 'LineWidth',LWidth, 'LineStyle',LStyle, 'Color',LColor, 'EraseMode',EraseMode );
    c = c + 1;
else
    h(c)    = line( [x-0 x+0],ylim, 'LineWidth',LWidth, 'LineStyle',LStyle, 'Color',LColor, 'EraseMode',EraseMode );
    c = c + 1;
    h(c)    = line( xlim,[y-0 y+0], 'LineWidth',LWidth, 'LineStyle',LStyle, 'Color',LColor, 'EraseMode',EraseMode );
    c = c + 1;
end
if ~strcmp(TLLabel,'')
    h(c) = text(x,y,[TLLabel ' '],...
        'EraseMode',EraseMode,...
        'HorizontalAlignment','right',...
        'VerticalAlignment','bottom',...
        'Color',TColor,...
        'FontSize',TSize,...
        'FontWeight',FontWeight);
    c   = c + 1;
end
if ~strcmp(TRLabel,'')
    h(c) = text(x,y,[' ' TRLabel],...
        'EraseMode',EraseMode,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','bottom',...
        'Color',TColor,...
        'FontSize',TSize,...
        'FontWeight',FontWeight);
    c   = c + 1;
end
if ~strcmp(BLLabel,'')
    h(c) = text(x,y+1,[BLLabel ' '],...
        'EraseMode',EraseMode,...
        'HorizontalAlignment','right',...
        'VerticalAlignment','top',...
        'Color',TColor,...
        'FontSize',TSize,...
        'FontWeight',FontWeight);
    c   = c + 1;
end
if ~strcmp(BRLabel,'')
    h(c) = text(x,y+1,[' ' BRLabel],...
        'EraseMode',EraseMode,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top',...
        'Color',TColor,...
        'FontSize',TSize,...
        'FontWeight',FontWeight);
    %c   = c + 1;
end
return;
end
%==========================================================================
% Drawing rectangle in current axes
%==========================================================================
function h = caDrawRectangle(X1,Y1,X2,Y2, EraseMode, LWidth,LColor,LStyle, FaceColor, TColor,TSize,FontWeight, TopLabel,BottomLabel,RightLabel,LeftLabel)
%==========================================================================
%   This functions draws rectangle with four labels in current axes
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% h = caDrawRectangle(X1,Y1,X2,Y2, EraseMode, LWidth,LColor,LStyle, FaceColor,...
%                       TColor,TSize,FontWeight, TopLabel,BottomLabel,RightLabel,LeftLabel)
% 
%       X1,X2,Y1,Y2     Coordinates of the rectangle
%       EraseMode       EraseMode property of rectangle ad all text labels
%       LWidth          Line width
%       LColor          Line color
%       LStyle          Line style
%       FaceColor       Face color
%       TColor          Text color
%       TSize           Text size
%       FontWeight      Font weight
%       TopLabel        Top label string
%       BottomLabel     Bottom label string
%       RightLabel      Right label string
%       LeftLabel       Left label string
% 
%       h               array of handles of all created elements
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Type
%           help rectangle
%   or
%           help text
%   for more information about rectangle and text parameters
%==========================================================================
c   = 1;
if (round(X1)~=round(X2))&&(round(Y1)~=round(Y2))
    h(c) = rectangle('Position',[X1,Y1,X2-X1,Y2-Y1],...
        'FaceColor',FaceColor,...
        'EraseMode',EraseMode,...
        'EdgeColor',LColor,...
        'LineWidth',LWidth,...
        'LineStyle',LStyle);
    c   = c + 1;
end
if ~strcmp(TopLabel,'')
    h(c) = text(X1,Y1,TopLabel,...
        'EraseMode',EraseMode,...
        'HorizontalAlignment','Left',...
        'VerticalAlignment','bottom',...
        'Color',TColor,...
        'FontSize',TSize,...
        'FontWeight',FontWeight);
    c   = c + 1;
end
if ~strcmp(BottomLabel,'')
    h(c) = text(X2,Y2+1,BottomLabel,...
        'EraseMode',EraseMode,...
        'HorizontalAlignment','Right',...
        'VerticalAlignment','top',...
        'Color',TColor,...
        'FontSize',TSize,...
        'FontWeight',FontWeight);
    c   = c + 1;
end
if ~strcmp(RightLabel,'')
    h(c) = text(X2,Y1,RightLabel,...
        'Rotation',-90.0,...
        'EraseMode',EraseMode,...
        'HorizontalAlignment','Left',...
        'VerticalAlignment','bottom',...
        'Color',TColor,...
        'FontSize',TSize,...
        'FontWeight',FontWeight);
    c   = c + 1;
end
if ~strcmp(LeftLabel,'')
    h(c) = text(X1,Y2,LeftLabel,...
        'Rotation',+90.0,...
        'EraseMode',EraseMode,...
        'HorizontalAlignment','Left',...
        'VerticalAlignment','bottom',...
        'Color',TColor,...
        'FontSize',TSize,...
        'FontWeight',FontWeight);
    %c   = c + 1;
end
return;

end
