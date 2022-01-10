function assemPotentialKandV(BEMesh,quadRuleA,quadRuleB)
	# Caution 
	# quadRuleA and B need to be different, otherwise during the double integration the Green function 
	# can be evaluated at the same points where it is singular


	#constant 
	mu0=4*pi*10^-7; mu0inv=1/mu0	

	# dimension
	n=size(BEMesh.edges,1)
	n2=n*n

	#matrix to be computed in sparse format
	VE=zeros(n2); VI=zeros(n2); VJ=zeros(n2) #P0-P0, one contribution by pair of edges

	KE1=zeros(n,2); KI1=zeros(n,2); KJ1=zeros(n,2) #P1-P0, two contributions by edge
	KE2=zeros(n2,2); KI2=zeros(n2,2); KJ2=zeros(n2,2) #P1-P0, two contributions by pair of edges
	#reshape at the end vector size(2*n+2*n2)

	#quadrature rules
	xGaussA=quadRuleA.x
	wGaussA=quadRuleA.w
	nQuadA=size(quadRuleA.w,1)
	xGaussB=quadRuleB.x
	wGaussB=quadRuleB.w
	nQuadB=size(quadRuleB.w,1)

	if(nQuadA==nQuadB)
		@show nQuadA
		@show nQuadB	
		println("The 2 quadratures list of points should be different for the double integration of the singular Green function")
		readline()
	end


	## single integral
	e=collect(1:n)
	edgeToNode=BEMesh.edges
	h=BEMesh.length
	for i in 1:nQuadA
		KE1=KE1+[0.5*h*0.5*(wGaussA[i]*(1-(0.5+xGaussA[i]*0.5)))  0.5*h*0.5*(wGaussA[i]*( (0.5+xGaussA[i]*0.5))) ]
		#The famous 1/2 for the angular coefficient is needed
		#However the BEM solution evaluated at a point not on the boundary but close to it is wrong, I guess because this anglar coefficient 
		# which should be 1 must be wrong if we are too close ...
		#without the 1/2:  
		#KE1=KE1+[h*0.5*(wGaussA[i]*(1-(0.5+xGaussA[i]*0.5)))  h*0.5*(wGaussA[i]*( (0.5+xGaussA[i]*0.5))) ]
	end

	KI1=[e e]
	KJ1=edgeToNode



	## double integrals
	# pairing of each edge with each other edge
	#ind=collect(1:n) #collect not needed
	ind=1:n
	e1=repeat(ind,1,n) 
	#e1=reshape(transpose(e1),n*n,1) # wrong, this gives a (n^2)x1 matrix
									# contrary to matlab a Nx1 matrix is not the same as a vector of size N !!! 
									# and then a dim is added in when using it as an index in edges
									# we need a vector of length n^2:
	e1=reshape(transpose(e1),n*n)

	#e2=repeat(ind,n,1)	# this gives a n^2x1 matrix and then a dim is added in when using it as index
						# we need a vector of length n^2:
	e2=repeat(ind,n)	
	#@show e1,e2

	#the two nodes of each edge
	edgeToNode1=BEMesh.edges[e1,:]
	edgeToNode2=BEMesh.edges[e2,:]
	#@show edgeToNode1
	#@show size(edgeToNode1)

	#in order to remap the quadrature we need the midpoints 
	#and length of first and second edge
    PA=0.5*(BEMesh.coordinates[edgeToNode1[:,1],:]+BEMesh.coordinates[edgeToNode1[:,2],:])	
    PB=0.5*(BEMesh.coordinates[edgeToNode2[:,1],:]+BEMesh.coordinates[edgeToNode2[:,2],:])
	#@show PA	

	hA=BEMesh.length[e1]
	hB=BEMesh.length[e2]
	#@show hA,hB

	#quadrature points on first and second edge
	PAG=zeros(size(PA,1),2,nQuadA)
	PBG=zeros(size(PB,1),2,nQuadB)
	for i in 1:nQuadA
		PAG[:,:,i]=PA + xGaussA[i]*0.5*(BEMesh.coordinates[edgeToNode1[:,2],:]-BEMesh.coordinates[edgeToNode1[:,1],:])
	end
	for i in 1:nQuadB
		PBG[:,:,i]=PB + xGaussB[i]*0.5*(BEMesh.coordinates[edgeToNode2[:,2],:]-BEMesh.coordinates[edgeToNode2[:,1],:])
	end

	#@show size(PAG)
	#@show PAG
	#@show size(PBG)
	#@show PBG

	#normal on the second edge
	normal2=BEMesh.normalsIn[e2,:]

	# pairing of each quadrature point of the first edge with each of the second edge
	idA=reshape(transpose(repeat(1:nQuadA,1,nQuadB)),nQuadA*nQuadB) 
	idB=repeat(1:nQuadB,nQuadA)	

	#@show idA
	#@show idB

	for i in 1:(nQuadA*nQuadB)
		#Evaluate Green function and first derivative (in y=B)
		G = @. Green(PAG[:,1,idA[i]],PAG[:,2,idA[i]],PBG[:,1,idB[i]],PBG[:,2,idB[i]],derivative=1)
		#G is a vector of named tuples, we need the broadcasted getindex function to access a value (with .() or @. see below)
		#@show G
		
		#Green function value
		G_P1GP2G = getindex.(G,:G_xy)
		#@show G_P1GP2G

		#Neumann trace
		dnG_P1GP2G = @.getindex(G,:dyrG_xy) .* normal2[:,1] +  @.getindex(G,:dyzG_xy) .* normal2[:,2]
		#@show dnG_P1GP2G 

		VE=VE+wGaussA[idA[i]]*wGaussB[idB[i]]*0.5*0.5*hA.*hB.*G_P1GP2G

		KE2=KE2+mu0inv*0.25*wGaussA[idA[i]]*wGaussB[idB[i]]*
			[ hA.*hB.*(1.0./PBG[:,1,idB[i]]).*dnG_P1GP2G*(1-(0.5+xGaussB[idB[i]]*0.5))  hA.*hB.*(1.0./PBG[:,1,idB[i]]).*dnG_P1GP2G*((0.5+xGaussB[idB[i]]*0.5)) ]

	end

	VI=e1; VJ=e2;
	KI2=[e1 e1]; KJ2=[edgeToNode2[:,1] edgeToNode2[:,2]]

	KI=[KI1;KI2]; KJ=[KJ1;KJ2]; KE=[KE1;KE2]
	KI=reshape(KI,size(KI,1)*size(KI,2))
	KJ=reshape(KJ,size(KJ,1)*size(KJ,2))
	KE=reshape(KE,size(KE,1)*size(KE,2))


	K=(E=KE,I=KI,J=KJ)
	V=(E=VE,I=VI,J=VJ)

	return K,V

end
