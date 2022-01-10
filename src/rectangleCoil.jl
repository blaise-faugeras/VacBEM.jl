function evalPsiFromRectangleCoil(r0,z0,dr,dz,Ic,quadRule_r,quadRule_z, xr,xz; derivative=0)
	
	surf=dr*dz
	jc=Ic/surf

	#define quadrature points and weigths in the rectangle
	xgr=r0 .+ quadRule_r.x *0.5 *dr	
	xgz=z0 .+ quadRule_z.x *0.5 *dz	

	xgr_rep=repeat(xgr,size(xgz,1))
	xgz_rep=reshape(repeat(xgz,1,size(xgr,1))',size(xgr,1)*size(xgz,1))
	
	wgr=quadRule_r.w *0.5 *dr	
	wgz=quadRule_z.w *0.5 *dz	

	wgr_rep=repeat(wgr,size(wgz,1))
	wgz_rep=reshape(repeat(wgz,1,size(wgr,1))',size(wgr,1)*size(wgz,1))

	#
	G = @.Green(xr,xz,xgr_rep,xgz_rep, derivative=derivative)
	G_PXPG = getindex.(G,:G_xy)

	val=wgr_rep.*wgz_rep.*G_PXPG*jc
	val=sum(val)	

	if(derivative==0)
		return (psi=val,)
	end
	if(derivative>0)
		dxrG = getindex.(G,:dxrG_xy)
		dxzG = getindex.(G,:dxzG_xy)

		dxrval=wgr_rep.*wgz_rep.*dxrG*jc
		dxzval=wgr_rep.*wgz_rep.*dxzG*jc
	
		dxrval=sum(dxrval)	
		dxzval=sum(dxzval)	
		
		return (psi=val, dxrpsi=dxrval, dxzpsi=dxzval)

	end

end

function assemblyRhsRectangleCoil(r0,z0,dr,dz,Ic,quadRule_r, quadRule_z, BEMesh, quadRule)

	# dimension
	n=size(BEMesh.edges,1)

	#
	rhs=zeros(n)

	#quadrature rule on BEMesh
	xGauss=quadRule.x
	wGauss=quadRule.w
	nQuad=size(quadRule.w,1)

	#in order to remap the quadrature we need the midpoints and length of edges
	e=collect(1:n)
	edgeToNode=BEMesh.edges
	h=BEMesh.length

    P=0.5*(BEMesh.coordinates[edgeToNode[:,1],:]+BEMesh.coordinates[edgeToNode[:,2],:])
	#@show P
	
	#quadrature points on all edges
	PG=zeros(size(P,1),2,nQuad)
	for i in 1:nQuad
		PG[:,:,i]=P + xGauss[i]*0.5*(BEMesh.coordinates[edgeToNode[:,2],:]-BEMesh.coordinates[edgeToNode[:,1],:])

		evalC=@.evalPsiFromRectangleCoil((r0,),(z0,),(dr,),(dz,),(Ic,),(quadRule_r,),(quadRule_z,), PG[:,1,i], PG[:,2,i])
		
		psiC = getindex.(evalC,:psi)

		rhs=rhs+0.5*wGauss[i]*h.*psiC
	end

	return rhs

end
