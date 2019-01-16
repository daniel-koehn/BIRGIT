function [grad] = median_filt(nx,ny,zeropad,zsmooth,grad)
%   median_filt Median filter
%   
%   Daniel Koehn
%   Kiel, the 23rd of June 2013
    
    % apply wavenumber domain damping
    % add boundary to the model to avoid filter artefacts on the inversion grid 
    ZIs = ones(ny+2.*zeropad,nx+2.*zeropad);

    for j=1:ny
    ZIs(j+zeropad,1:zeropad) = grad(j,1);
    ZIs(j+zeropad,nx+zeropad:nx+2.*zeropad) = grad(j,1);
    end

    ZIs(zsmooth+zeropad:ny+zeropad,1+zeropad:nx+zeropad) = grad(zsmooth:ny,1:nx); 

    for i=1:nx+(2.*zeropad)
    ZIs(1:zeropad+zsmooth,i) = ZIs(zsmooth+zeropad,i);
    ZIs(ny+zeropad:ny+2.*zeropad,i) = ZIs(ny+zeropad,i);
    end
    
    % apply median filter
    ZIs=medfilt2(ZIs,[zeropad zeropad]);
    
    grad(1:ny,1:nx)=ZIs(1+zeropad:ny+zeropad,1+zeropad:nx+zeropad);
    clear ZIs;

end

