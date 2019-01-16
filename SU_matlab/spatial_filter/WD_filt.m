function [grad] = WD_filt(nx,ny,zeropad,zsmooth,grad,damp)
%   WD_filt Wavenumber domain filter
%   
%   Daniel Koehn
%   Kiel, the 10th of April 2013
    
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
    
    nxtmp = nx + 2.*zeropad;
    nytmp = ny + 2.*zeropad;
    
    kxc = round(nxtmp./2)+1;
    kzc = round(nytmp./2)+1;

    fftvp3 = fft2(ZIs);
    fftvp3 = fftshift(fftvp3);

    kx = 1:nxtmp;
    kz = 1:nytmp;

    [KX,KZ] = meshgrid(1:nxtmp,1:nytmp);

    gausstaper = exp(-damp.*((KX-kxc).^2+(KZ-kzc).^2));
    fftvp3 = fftvp3.*gausstaper;

    ifftvp3 = ifftshift(fftvp3);
    ifftvp3 = ifft2(ifftvp3);
    ZIs = real(ifftvp3);
    
    grad(1:ny,1:nx)=ZIs(1+zeropad:ny+zeropad,1+zeropad:nx+zeropad);
    clear ZIs;

end

