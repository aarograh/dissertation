classdef quadratureClass
    %QUADRATURECLASS Contains quadrature information
    %   Contains the angles, their cosines, and weights
    %   for a polar quadrature for 1D MOC
    
    properties
        npol
        angles
        cosines
        weights
    end
    
    methods
        function obj = quadratureClass(npol)
            %QUADRATURECLASS Sets up a quadrature class
            %   npol - Number of polar angles.  Currently, only values of
            %          1 are accepted.
            
            obj.npol = npol/2;
            [obj.cosines, obj.weights] = createQuadrature(npol);
            obj.angles = acos(obj.cosines);
        end
    end
    
end

function [mu w] = createQuadrature(N)
   %N is ordinate set: 2, 4, 8, ...
   %Return symmetric gauss-legendre quadrature set
   % syntax: [mu w] = createQuadrature( 4 )
   if (N == 2)
      mu = [.577350269189626]; %#ok<NBRAK>
      w  = [1.0]; %#ok<NBRAK>
   elseif (N == 4)
      mu = [.339981043584856 .861136311594053];
      w  = [.652145154862546 .347854845137454];
   elseif (N == 8)
      mu = [.183434642495650 .525532409916329 .796666477413627 .960289856497536];
      w  = [.362683783378363 .313706645877887 .222381034453374 .101228536290376];
   elseif (N == 16)
      mu = [0.989400934991650,0.944575023073233,0.865631202387832,...
         0.755404408355003,0.617876244402644,0.458016777657227,...
         0.281603550779259,0.0950125098376370;];
      w = [0.0271524594117540,0.0622535239386480,0.0951585116824930,...
         0.124628971255534,0.149595988816577,0.169156519395003,...
         0.182603415044924,0.189450610455067;];
   elseif (N == 32)
      mu = [0.0483076656877380,0.144471961582796,0.239287362252137,...
         0.331868602282128,0.421351276130635,0.506899908932229,...
         0.587715757240762,0.663044266930215,0.732182118740290,...
         0.794483795967942,0.849367613732570,0.896321155766052,...
         0.934906075937740,0.964762255587506,0.985611511545268,...
         0.997263861849482];
      w = [0.0965400885147260,...
         0.0956387200792750,0.0938443990808050,0.0911738786957640,...
         0.0876520930044040,0.0833119242269470,0.0781938957870700,...
         0.0723457941088490,0.0658222227763620,0.0586840934785360,...
         0.0509980592623760,0.0428358980222270,0.0342738629130210,...
         0.0253920653092620,0.0162743947309060,0.007018610009470];
   else
      error('Sorry, invalid N for quadrature')
   end
   
   assert(length(mu) == length(w));
   assert(abs(sum(w) - 1.0) < 2*eps);
   assert(all(mu > 0.0) && all(mu < 1.0));
   assert(all(w > 0));
end