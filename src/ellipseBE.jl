function ellipseBE(x0,y0,a,b,n)

	#==
	n points defining the ellipse ((x-x0)/a)^2+((y-y0)/b)^2=1
	counter-clockwise
	==#
	dtheta=2pi/n
	coordinates=zeros(n,2)
	for i in 1:n
		coordinates[i,1]=x0+a*cos((i-1)*dtheta)
		coordinates[i,2]=y0+b*sin((i-1)*dtheta)
	end

	# corresponding edges
	edges=[i+j for i in 1:n, j in 0:1]
	edges[n,2]=1

	# and length and inner normals
	length=zeros(n)
	normalsIn=zeros(n,2)
	for i in 1:n
		tx=coordinates[edges[i,2],1]-coordinates[edges[i,1],1]
		ty=coordinates[edges[i,2],2]-coordinates[edges[i,1],2]
		length[i]=norm([tx ty])
		normalsIn[i,:]=[-ty tx]/norm([tx ty])
	end
	
	

	return (coordinates=coordinates, edges=edges, length=length, normalsIn=normalsIn)

end
