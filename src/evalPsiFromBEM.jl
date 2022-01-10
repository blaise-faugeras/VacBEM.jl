function evalPsiFromBEM(BEMesh,quadRule,f,q,xr,xz; derivative=0)
#=
Given a BE contour BEMesh, a quadrature rule on this contour, the Dirichlet P1 dof and the Neumann P0 dof, 
evaluate psi at any point (xr,xz) outside the contour resulting from BEM
q is a vector of size n edges
f is a vector of size n nodes (the BEMesh contour is closed)
=#

	#constant 
	mu0=4*pi*10^-7; mu0inv=1/mu0	
	
	# dimension
	n=size(BEMesh.edges,1)

	#quadrature rule
	xGauss=quadRule.x
	wGauss=quadRule.w
	nQuad=size(quadRule.w,1)

	#in order to remap the quadrature we need the midpoints and length of edges
	e=collect(1:n)
	edgeToNode=BEMesh.edges
	h=BEMesh.length

    P=0.5*(BEMesh.coordinates[edgeToNode[:,1],:]+BEMesh.coordinates[edgeToNode[:,2],:])
	#@show P
	
	#normal 
	normal=BEMesh.normalsIn
	
	#quadrature points on all edges
	PG=zeros(size(P,1),2,nQuad)
	for i in 1:nQuad
		PG[:,:,i]=P + xGauss[i]*0.5*(BEMesh.coordinates[edgeToNode[:,2],:]-BEMesh.coordinates[edgeToNode[:,1],:])
	end

	#evaluate Green function at all y=Gauss points PG and x=(xr,xz)
	#need to broadcast on the number of edges
	
	DE=zeros(n,2); DI=zeros(n,2); DJ=zeros(n,2)
	NE=zeros(n); NI=zeros(n); NJ=zeros(n)
	
	dxrDE=zeros(n,2); 
	dxzDE=zeros(n,2); 
	dxrNE=zeros(n)
	dxzNE=zeros(n)

	for i in 1:nQuad
		#Evaluate Green function and first derivative (in y=B)	
		G = @. Green(xr,xz,PG[:,1,i],PG[:,2,i],derivative=derivative+1)  
		#G is a vector of named tuples, we need the broadcasted getindex function to access a value (with .() or @. see below)
		#@show G
		
		#Green function value
		G_PXPG = getindex.(G,:G_xy)
		#@show G_PXPG

		#Neumann trace
		dnG_PXPG = @.getindex(G,:dyrG_xy) .* normal[:,1] +  @.getindex(G,:dyzG_xy) .* normal[:,2]
		#@show dnG_P1GP2G 

		# Matrices. Matrix.vector product to be summed afterwards for integral values
		# for psi(x)
		NE=NE+0.5*wGauss[i]*h.*G_PXPG
	
		#DE=DE-0.5*wGauss[i]*[ h.*dnG_PXPG*(1-(0.5+xGauss[i]*0.5))  h.*dnG_PXPG*((0.5+xGauss[i]*0.5)) ]
		DE=DE-0.5*wGauss[i]*mu0inv*[ h.*(1.0./PG[:,1,i]).*dnG_PXPG*(1-(0.5+xGauss[i]*0.5))  h.*(1.0./PG[:,1,i]).*dnG_PXPG*((0.5+xGauss[i]*0.5)) ]

		# for d_xr_psi(x) and d_xz_psi(x)
		if(derivative>0)
			dxrG= @.getindex(G,:dxrG_xy);
			dxzG= @.getindex(G,:dxzG_xy);
			dxrNE=dxrNE+0.5*wGauss[i]*h.*dxrG
			dxzNE=dxzNE+0.5*wGauss[i]*h.*dxzG

			dxrdyn = @.getindex(G,:dxrdyrG_xy) .* normal[:,1] + @.getindex(G,:dxrdyzG_xy) .* normal[:,2]
			dxzdyn = @.getindex(G,:dxzdyrG_xy) .* normal[:,1] + @.getindex(G,:dxzdyzG_xy) .* normal[:,2]

			dxrDE=dxrDE-0.5*wGauss[i]*mu0inv*[ h.*(1.0./PG[:,1,i]).*dxrdyn*(1-(0.5+xGauss[i]*0.5))  h.*(1.0./PG[:,1,i]).*dxrdyn*((0.5+xGauss[i]*0.5)) ]
			dxzDE=dxzDE-0.5*wGauss[i]*mu0inv*[ h.*(1.0./PG[:,1,i]).*dxzdyn*(1-(0.5+xGauss[i]*0.5))  h.*(1.0./PG[:,1,i]).*dxzdyn*((0.5+xGauss[i]*0.5)) ]
		end




	end

	NI=e; NJ=e
	DI=[e e]; DJ=edgeToNode

	DI=reshape(DI,size(DI,1)*size(DI,2))
	DJ=reshape(DJ,size(DJ,1)*size(DJ,2))
	DE=reshape(DE,size(DE,1)*size(DE,2))

	Dspmat=sparse(DI,DJ,DE,n,n)
	Nspmat=sparse(NI,NJ,NE,n,n)
	Ival=Dspmat*f+Nspmat*q
	Ival=sum(Ival)

	if(derivative==0)
		return (psi=Ival,)
	end

	if(derivative>0)
		dxrDE=reshape(dxrDE,size(dxrDE,1)*size(dxrDE,2))
		dxzDE=reshape(dxzDE,size(dxzDE,1)*size(dxzDE,2))
		dxrDspmat=sparse(DI,DJ,dxrDE,n,n)
		dxzDspmat=sparse(DI,DJ,dxzDE,n,n)
		dxrNspmat=sparse(NI,NJ,dxrNE,n,n)
		dxzNspmat=sparse(NI,NJ,dxzNE,n,n)
		
		dxrIval=dxrDspmat*f+dxrNspmat*q
		dxrIval=sum(dxrIval)
		
		dxzIval=dxzDspmat*f+dxzNspmat*q
		dxzIval=sum(dxzIval)

		return (psi=Ival, drpsi=dxrIval, dzpsi=dxzIval)
	end

	
end
