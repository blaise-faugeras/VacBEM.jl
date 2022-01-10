push!(LOAD_PATH,pwd()*"/src") #using local modules
using Revise #enables to change the source code and not having to reload all modules or kill the REPL
using Plots
using SparseArrays

using QuadratureRules
#using GreenDeltaStar
using MeshFactory
using Assembly

#ellipse
x0=2.0; y0=0.0; a=0.1; b=0.5
n=30
BEMesh=ellipseBE(x0,y0,a,b,n)

#p=plot(BEMesh.coordinates[:,1],BEMesh.coordinates[:,2], aspect_ratio=:equal)
#display(p)

#quadrature rules
quadRuleA=gauleg(-1,1,3)
#quadRuleA=gauleg(-1,1,5)
quadRuleB=gauleg(-1,1,5)

#assemble BEM matrices
#@time K,V=assemPotentialKandV(BEMesh,quadRuleA,quadRuleB)
@time K,V=assemPotentialKandVwithDuffy(BEMesh,quadRuleA,quadRuleB)
#@show K
#@show V

#sparse matrix
Kspmat=sparse(K.I,K.J,K.E,n,n)
#@show Kspmat

Vspmat=sparse(V.I,V.J,V.E) #no need to set the dimension of the matrix, sparse seeks the largest i and j
#@show Vspmat

#coils 
quadRule_r=gauleg(-1,1,4)
quadRule_z=gauleg(-1,1,4)
quadRule=gauleg(-1,1,2)

r0=[3,3,0.5]; z0=[3,-3,0]; dr=[0.25,0.25,0.2]; dz=[0.25,0.25,2]
Ic=[5.0e5, 5.0e5,-5.0e6]
#Ic=zeros(3)

Jcrhs=zeros(n)

for i=1:size(Ic,1)
	rhs=assemblyRhsRectangleCoil(r0[i],z0[i],dr[i],dz[i],Ic[i],quadRule_r, quadRule_z, BEMesh, quadRule)
	global Jcrhs=Jcrhs+rhs
end

@show Jcrhs



##boundary integral equation K.f-V.q = Jc with f P1 Dirichlet and q P0 Neumann
f=2.0*ones(n) # Dirichlet given
@time q=Vspmat\(Kspmat*f - Jcrhs) #compute Neumann
#boundary integral equation K.f-V.q=0 with f P1 Dirichlet and q P0 Neumann
#q=ones(n)*1.0e6 # Neumann given
#f=Kspmat\(Vspmat*q) #compute Dirichlet

@show f
@show q

##
##
##xr=[2,2.02,1.92]; xz=[0,0,0]
###psi=evalPsiFromBEM(BEMesh,quadRule,f,q,xr,xz)
##psi=@.evalPsiFromBEM((BEMesh,),(quadRule,),(f,),(q,),xr,xz) 
###@show psi

ln=50
xr=range(0.1,6,length=ln)
#xr=range(1.5,2.5,length=ln)
xr_rep=repeat(xr,ln)

xz=range(-5,5,length=ln)
#xz=range(-1.5,1.5,length=ln)
xz_rep=reshape(repeat(xz,1,ln)',ln*ln)


@time evalBEM=@.evalPsiFromBEM((BEMesh,),(quadRule,),(f,),(q,),xr_rep,xz_rep) 
psiBEM=@.getindex(evalBEM,:psi)
#@show psiBEM



# points inside ellipse are included in plot but psi is wrong of course 
#pBEM=contourf(xr,xz,psiBEM)
#display(pBEM)


## coils
#quadRule_r=gauleg(-1,1,4)
#quadRule_z=gauleg(-1,1,4)
#
#r0=[3,3,0.5]; z0=[3,-3,0]; dr=[0.25,0.25,0.2]; dz=[0.25,0.25,2]
#Ic=[5.0e5, 5.0e5,-5.0e6]

#xr=4; xz=4;
#psiC=evalPsiFromRectangleCoil(r0,z0,dr,dz,Ic,quadRule_r,quadRule_z,xr,xz)

psiC=zeros(size(xr_rep,1))

for i=1:size(Ic,1)
	local evalC=@.evalPsiFromRectangleCoil((r0[i],),(z0[i],),(dr[i],),(dz[i],),(Ic[i],),(quadRule_r,),(quadRule_z,),xr_rep,xz_rep)
	global psiC=psiC+@.getindex(evalC,:psi)
end
#@show psiC

#pC=contourf(xr,xz,psiC)
#display(pC)


psi=psiBEM+psiC
@time p=contourf(xr,xz,psi)
#p=contour(xr,xz,psi)
display(p)

dl=0.01
xr=[x0-a-dl,x0+a+dl,x0,x0]; xz=[y0,y0,y0+b+dl,y0-b-dl]
evalBEM=@.evalPsiFromBEM((BEMesh,),(quadRule,),(f,),(q,),xr,xz,derivative=0,) #trick to broadcast, "fixed" arguments are seen as tuples so almost a scalar
psiBEM=@.getindex(evalBEM,:psi)

psiC=zeros(size(xr,1))
for i=1:size(Ic,1)
	local evalC=@.evalPsiFromRectangleCoil((r0[i],),(z0[i],),(dr[i],),(dz[i],),(Ic[i],),(quadRule_r,),(quadRule_z,),xr,xz)
	global psiC=psiC+@.getindex(evalC,:psi)
end
psi=psiBEM+psiC
@show psi


println("This is the end")
readline()



