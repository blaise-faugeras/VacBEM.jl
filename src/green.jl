function Green(xr,xz,yr,yz; derivative=0)
    #= 
    Fundamental solution (Green function) and derivatives of the Grad-Shafranov Delta^* operator
    -Delta_x^* G(x,y)=delta(x-y)
    y=(yr,yz) is the source point (a filament with unit current)
    x=(xr,xz) is the evaluation point
    Derivatives from Itagaki, Engineering Analysis with Boundary Elements, 33 (2009) 845-857
    
    derivative=0 only Green function value
    derivative=1 + gradient wrt to x and y
    derivative=2 + second order cross derivatives (see return outputs below)
    =#
    
    mu0=4*pi*10^-7
    #modulus at points x and y
    k_xy=sqrt(4*xr*yr/((xr+yr)^2+(xz-yz)^2))
    
    #Complete elliptic integrals of first and second kind
    K_xy,E_xy=Elliptic.ellipke(k_xy^2)
    
    #Green function
    G_xy = (mu0/pi*((xr*yr)^(1/2))/k_xy)*((1-0.5*k_xy^2)*K_xy - E_xy)
    
    
    #println("derivative=",derivative)
    
    if(derivative==0)
        return (G_xy=G_xy,)
        #return G_xy
    end
       
    if(derivative>0)
        # gradient with respect to y
        dyrG_xy = (mu0/(2*pi))*((yr   )/(sqrt( (xr+yr)^2 + (xz-yz)^2 )))*(K_xy - ((yr^2-xr^2+(xz-yz)^2 )/( (xr-yr)^2 + (xz-yz)^2))*E_xy)
        dyzG_xy = (mu0/(2*pi))*((yz-xz)/(sqrt( (xr+yr)^2 + (xz-yz)^2 )))*(K_xy - ((yr^2+xr^2+(xz-yz)^2 )/( (xr-yr)^2 + (xz-yz)^2))*E_xy)
    
        # gradient with respect to x
        dxrG_xy = (mu0/(2*pi))*((xr   )/(sqrt( (xr+yr)^2 + (xz-yz)^2 )))*(K_xy - ((xr^2-yr^2+(xz-yz)^2 )/( (xr-yr)^2 + (xz-yz)^2))*E_xy)
        dxzG_xy = (mu0/(2*pi))*((xz-yz)/(sqrt( (xr+yr)^2 + (xz-yz)^2 )))*(K_xy - ((xr^2+yr^2+(xz-yz)^2 )/( (xr-yr)^2 + (xz-yz)^2))*E_xy)
        
        if(derivative==1)
            #return G_xy, dyrG_xy, dyzG_xy, dxrG_xy, dxzG_xy
            return (G_xy=G_xy, dyrG_xy=dyrG_xy, dyzG_xy=dyzG_xy, dxrG_xy=dxrG_xy, dxzG_xy=dxzG_xy)
        end
        
    end
    
    if(derivative==2)
    
        # second order cross derivatives
        dxrdyrG_xy = -((xr+yr)/((xr+yr)^2+(xz-yz)^2))*dxrG_xy + 
            (mu0/2*pi)*xr*(1/((xr+yr)^2+(xz-yz)^2)^(1/2))*(1/((xr-yr)^2+(xz-yz)^2))*
            ( ((xr^2-yr^2+(xz-yz)^2)/((xr+yr)^2+(xz-yz)^2))*((xr-yr)*K_xy+(xr+yr)*E_xy) -2*xr*E_xy + ((4*yr*(xz-yz)^2)/((xr-yr)^2+(xz-yz)^2))*E_xy )
    
        dxzdyrG_xy = -((xr+yr)/((xr+yr)^2+(xz-yz)^2))*dxzG_xy + 
            (mu0/2*pi)*xr*(xz-yz)*(1/((xr-yr)^2+(xz-yz)^2)^(1/2))*(1/((xr-yr)^2+(xz-yz)^2))*
            ( ((xr^2-yr^2+(xz-yz)^2)/((xr+yr)^2+(xz-yz)^2))*(K_xy+E_xy) -2*E_xy - ((4*yr*(xr-yr))/((xr-yr)^2+(xz-yz)^2))*E_xy )
    
        dxrdyzG_xy = ((xz-yz)/((xr+yr)^2+(xz-yz)^2))*dxrG_xy +
            (mu0/pi)*yr*xr*(xz-yz)*
            (1/((xr+yr)^2+(xz-yz)^2)^(1/2))*(1/((xr-yr)^2+(xz-yz)^2))*
            ( (1/((xr+yr)^2+(xz-yz)^2))*((xr-yr)*K_xy+(xr+yr)*E_xy) +
            ((-2*(xr-yr))/((xr-yr)^2+(xz-yz)^2))*E_xy )
    
        dxzdyzG_xy = (((-1)/(xz-yz))+((xz-yz)/((xr+yr)^2+(xz-yz)^2)))*dxzG_xy +
            (mu0/pi)*yr*xr*(xz-yz)^2 *
            (1/((xr+yr)^2+(xz-yz)^2)^(1/2))*(1/((xr-yr)^2+(xz-yz)^2))*
            ( (1/((xr+yr)^2+(xz-yz)^2))*(K_xy+E_xy) +
            (-2/((xr-yr)^2+(xz-yz)^2))*E_xy )
    
        #return G_xy, dyrG_xy, dyzG_xy, dxrG_xy, dxzG_xy, dxrdyrG_xy, dxzdyrG_xy, dxrdyzG_xy, dxzdyzG_xy
        return (G_xy=G_xy, dyrG_xy=dyrG_xy, dyzG_xy=dyzG_xy, dxrG_xy=dxrG_xy, dxzG_xy=dxzG_xy, dxrdyrG_xy=dxrdyrG_xy, dxzdyrG_xy=dxzdyrG_xy, dxrdyzG_xy=dxrdyzG_xy, dxzdyzG_xy=dxzdyzG_xy)
    end

end
