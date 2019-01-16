function [grad] = gauss_filt(nx,ny,zeropadx,zeropady,zsmooth,grad)
%   gauss_filt Gauss filter
%   
%   Daniel Koehn
%   Kiel, the 23rd of June 2013
    
    % apply wavenumber domain damping
    % add boundary to the model to avoid filter artefacts on the inversion grid 
    ZIs = ones(ny+2.*zeropady,nx+2.*zeropadx);

    grad_old = grad;
    
    for j=1:ny
    ZIs(j+zeropady,1:zeropadx) = grad(j,1);
    ZIs(j+zeropady,nx+zeropadx:nx+2.*zeropadx) = grad(j,1);
    end

    ZIs(zsmooth+zeropady:ny+zeropady,1+zeropadx:nx+zeropadx) = grad(zsmooth:ny,1:nx); 

    for i=1:nx+(2.*zeropadx)
    ZIs(1:zeropady+zsmooth,i) = ZIs(zsmooth+zeropady,i);
    ZIs(ny+zeropady:ny+2.*zeropady,i) = ZIs(ny+zeropady,i);
    end
    
    % apply gaussian filter
    H=fspecial('gaussian',zeropady,zeropadx);
    ZIs = imfilter(ZIs,H,'replicate');
    
    grad(1:ny,1:nx)=ZIs(1+zeropady:ny+zeropady,1+zeropadx:nx+zeropadx);
    grad(1:zsmooth,1:nx)=grad_old(1:zsmooth,1:nx);
    
    clear ZIs;

end

