import Pkg; Pkg.add("ITensors")

import Pkg; Pkg.add("Plots")

import Pkg; Pkg.add("ITensorMPS")

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  tau = 0.1
  ttotal = 200
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HEgamma0.4_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 4.0 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HEgamma4.0_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=255
  tau = 0.1
  ttotal = 200
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HEgamma0.4_128.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=255
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=4.0
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HEgamma4.0_128.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 511
  tau = 0.1
  ttotal = 200
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HEgamma0.4_256.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 511
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 4.0
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HEgamma4.0_256.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  tau = 0.1
  ttotal = 200
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = (m_0+h)*op("Sz",s1)
     p=1
    elseif p==1
     hj = -(m_0+h)*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_numopgamma0.4_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  tau = 0.1
  ttotal = 200
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=4.0
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = (m_0+h)*op("Sz",s1)
     p=1
    elseif p==1
     hj = -(m_0+h)*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_numopgamma4.0_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 255
  tau = 0.1
  ttotal = 200
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = (m_0+h)*op("Sz",s1)
     p=1
    elseif p==1
     hj = -(m_0+h)*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_numopgamma0.4_128.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 255
  tau = 0.1
  ttotal = 200
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 4.0
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = (m_0+h)*op("Sz",s1)
     p=1
    elseif p==1
     hj = -(m_0+h)*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_numopgamma4.0_128.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 511
  tau = 0.1
  ttotal = 200
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = (m_0+h)*op("Sz",s1)
     p=1
    elseif p==1
     hj = -(m_0+h)*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_numopgamma0.4_256.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 511
  tau = 0.1
  ttotal = 200
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 4.0 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = (m_0+h)*op("Sz",s1)
     p=1
    elseif p==1
     hj = -(m_0+h)*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_numopgamma4.0_256.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

N=15
b=(N/2)
for j in 1:2:N-2
 if 1<=j<= b-2
    println(j)
 end
end

for j in 1:2:N-2
 if b-2<=j<= N-1
    println(j)
 end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  b=(N/2)
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.5
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
    if 1<=j<= b-2    
     hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
     Gj = exp(-im * (tau) *hj)
     push!(gates1, Gj)
    elseif b-2<j<= N-2
      hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
      Gj = exp(-im * (tau) *hj)
     push!(gates1, Gj)
    end
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  bond_dim_max=1000
  filename = "SvN_HIgamma0.5_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  b=(N/2)
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIfullgamma0.4_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=255
  b=(N/2)
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIfullgamma0.4_128.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=511
  b=(N/2)
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIfullgamma0.4_256.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  b=(N/2)
  tau = 0.1
  ttotal = 80
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 1.5 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIfullgamma1.5_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f= 1
  N= 255
  b=(N/2)
  tau = 0.1
  ttotal = 80
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 1.5 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIfullgamma1.5_128.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f= 1
  N= 511
  b=(N/2)
  tau = 0.1
  ttotal = 80
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 1.5 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIfullgamma1.5_256.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

N=15
b=(N/2)
for j in 1:2:N-2
    if (1 <= j <= b-2) || (b<= j <= N-2)
        println(j)
    end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f= 1
  N= 127
  b=(N/2)
  tau = 0.1
  ttotal = 80
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 0.5 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
    if (1 <= j <= b-2) || (b<= j <= N-2)
     hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    else
      hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    end
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIlinkexgamma0.5_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f= 1
  N= 127
  b=(N/2)
  tau = 0.01
  ttotal = 20
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 0.5 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
    if (1 <= j <= b-2) || (b<= j <= N-2)
     hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    else
      hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    end
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIlinkexgamma0.5_64_tau0.01.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f= 1
  N= 255
  b=(N/2)
  tau = 0.1
  ttotal = 80
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 0.5 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
    if (1 <= j <= b-2) || (b<= j <= N-2)
     hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    else
      hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    end
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIlinkexgamma0.5_128.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f= 1
  N= 511
  b=(N/2)
  tau = 0.1
  ttotal = 80
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 0.5 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
    if (1 <= j <= b-2) || (b<= j <= N-2)
     hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    else
      hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    end
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIlinkexgamma0.5_256.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f= 1
  N= 127
  b=(N/2)
  tau = 0.1
  ttotal = 80
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 4.0 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
    if (1 <= j <= b-2) || (b<= j <= N-2)
     hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    else
      hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    end
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIlinkexgamma4.0_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f= 1
  N= 255
  b=(N/2)
  tau = 0.1
  ttotal = 80
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 4.0 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
    if (1 <= j <= b-2) || (b<= j <= N-2)
     hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    else
      hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    end
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIlinkexgamma4.0_128.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f= 1
  N= 511
  b=(N/2)
  tau = 0.1
  ttotal = 80
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 4.0 
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
    if (1 <= j <= b-2) || (b<= j <= N-2)
     hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    else
      hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    end
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIlinkexgamma4.0_256.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  b=(N/2)
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 1.5
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
    if 1<=j<= b-2
     hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
     Gj = exp(-im * (tau) *hj)
     push!(gates1, Gj)
    elseif b-2<j<= N-2
      hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
      Gj = exp(-im * (tau / 2) *hj)
      push!(gates1, Gj)
    end
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIsubgamma1.5_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  b=(N/2)
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 2.0
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
    if 1<=j<= b-2
     hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
     Gj = exp(-im * (tau) *hj)
     push!(gates1, Gj)
    elseif b-2<j<= N-1
      hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
      Gj = exp(-im * (tau / 2) *hj)
      push!(gates1, Gj)
    end
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIsubgamma2.0_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=15
  cutoff = 1E-10
  tau = 0.1
  ttotal = 30
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates=append!(gates1,gates3,gates4)
  #gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  nsweeps = 5
  maxdim = [500,1200]
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename ="EEnh_check_TEBD1storder.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=15
  tau = 0.1
  ttotal = 30
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "EEnh_check_TEBD2ndorder.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

let
    ############################
    # PARAMETERS
    ############################
    nsites  = 8
    N       = 2*nsites - 1
    b       = N ÷ 2
    tau     = 0.1
    ttotal  = 30.0
    X       = 0.5
    gamma   = -im*0.4
    f       = 1.0
    mu      = 2*sqrt(X)
    cutoff  = 1e-8
    maxdim  = 200

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    p = 0
    state = String[]

    for n in 1:N
        if iseven(n)
            push!(state, "Dn")
        else
            if p == 0
                push!(state, "Dn")
                p = 1
            else
                push!(state, "Up")
                p = 0
            end
        end
    end

    psi0 = MPS(s, state)

    ############################
    # BUILD HAMILTONIAN MPO
    ############################
    ampo = AutoMPO()

    # Gauge-invariant hopping
    for j in 1:2:N-2
        add!(ampo, X, "S+", j, "Sx", j+1, "S-", j+2)
        add!(ampo, X, "S-", j, "Sx", j+1, "S+", j+2)
    end

    # Electric field
    for j in 2:2:N-1
        add!(ampo, (f + gamma), "Sz", j)
    end

    # Mass term
    for j in 1:2:N
        k = (j - 1) ÷ 2
        μ = (-1)^k * mu
        add!(ampo, μ, "Sz", j)
    end

    H = MPO(ampo, s)

    ############################
    # TDVP TIME EVOLUTION
    ############################
    println("# t    <E>")
    filename = "EEnh_check_TDVP.dat"

    open(filename, "w") do io
        for t in 0:tau:ttotal
            psi0 = psi0 + 1e-6 * randomMPS(s; linkdims=2)
            S_vN= ent1(psi0)
            #ExpE_t = exp_electric_field(psi0, s; f=1.0)
            println("$t $S_vN")
            psi0 = tdvp(H, -im * tau, psi0;
                        nsweeps=2, maxdim=maxdim,
                        cutoff=cutoff, nsite=1,
                        normalize=true)
            write(io, "$t \t $S_vN\n")
        end
    end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
  ############################
  # PARAMETERS
  ############################
  f = 1
  N = 15
  tau = 0.1
  ttotal = 10
  cutoff = 1E-8
  bond_dim_max = 200

  s = siteinds("S=1/2", N)

  # Initial staggered state
  state = String[]
  alt = ["Dn", "Up"]
  alt_index = 1
  for i in 1:N
    if iseven(i)
      push!(state, "Dn")
    else
      push!(state, alt[alt_index])
      alt_index = alt_index % 2 + 1
    end
  end

  x = 0.5
  m_0 = 2 * sqrt(x)
  #h_1 = 0.4

  gamma_list = [0.1, 0.4]   # <<< two gamma value

  for gamma in gamma_list

    println("\nRunning for gamma = $gamma")

    h = -im * gamma
    psi0 = MPS(s, state)   # reset state for each gamma

    gates1 = ITensor[]
    gates3 = ITensor[]
    gates4 = ITensor[]

    ################# Hopping term #################
    for j in 1:2:N-2
      s1 = s[j]
      s2 = s[j+2]
      s3 = s[j+1]

      hj = x * (
        op("S+", s1) * op("Sx", s3) * op("S-", s2) +
        op("S+", s2) * op("Sx", s3) * op("S-", s1)
      )

      Gj = exp(-im * tau * hj)
      push!(gates1, Gj)
    end

    ################# Electric field #################
    for j in 2:2:N-1
      hj = f * op("Sz", s[j])
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates3, Gj)
    end

    ################# Mass + non-Hermitian #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0 + h) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0 + h) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)
    end

    ################# Full Trotter step #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# Output #################
    filename = "SvN_check_$(gamma).dat"
    open(filename, "w") do io
      for t in 0:tau:ttotal
        S_vN = ent1(psi0)
        write(io, "$t\t$S_vN\n")
        println("gamma=$gamma  t=$t  SvN=$S_vN")

        psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
        normalize!(psi0)
      end
    end
  end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
  ############################
  # PARAMETERS
  ############################
  f = 1
  N = 127
  tau = 0.1
  ttotal = 100
  cutoff = 1E-8
  bond_dim_max = 1000

  s = siteinds("S=1/2", N)

  # Initial staggered state
  state = String[]
  alt = ["Dn", "Up"]
  alt_index = 1
  for i in 1:N
    if iseven(i)
      push!(state, "Dn")
    else
      push!(state, alt[alt_index])
      alt_index = alt_index % 2 + 1
    end
  end

  x = 0.5
  m_0 = 2 * sqrt(x)
  #h_1 = 0.4

  gamma_list1 = [0.1, 0.5,1]   # <<< two gamma value

  for gamma in gamma_list1

    println("\nRunning for gamma = $gamma")

    h = -im * gamma
    psi0 = MPS(s, state)   # reset state for each gamma

    gates1 = ITensor[]
    gates3 = ITensor[]
    gates4 = ITensor[]

    ################# Hopping term #################
    for j in 1:2:N-2
      s1 = s[j]
      s2 = s[j+2]
      s3 = s[j+1]

      hj = (x+h) * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

      Gj = exp(-im * tau * hj)
      push!(gates1, Gj)
    end

    ################# Electric field #################
    for j in 2:2:N-1
      hj = f * op("Sz", s[j])
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates3, Gj)
    end

    ################# Mass + non-Hermitian #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)
    end

    ################# Full Trotter step #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# Output #################
    filename = "SvN_HIfullgamma$(gamma)_64.dat"
    open(filename, "w") do io
      for t in 0:tau:ttotal
        S_vN = ent1(psi0)
        write(io, "$t\t$S_vN\n")
        println("gamma=$gamma  t=$t  SvN=$S_vN")

        psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
        normalize!(psi0)
      end
    end
  end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
  ############################
  # PARAMETERS
  ############################
  f = 1
  N = 127
  tau = 0.1
  ttotal = 100
  cutoff = 1E-8
  bond_dim_max = 1000

  s = siteinds("S=1/2", N)

  # Initial staggered state
  state = String[]
  alt = ["Dn", "Up"]
  alt_index = 1
  for i in 1:N
    if iseven(i)
      push!(state, "Dn")
    else
      push!(state, alt[alt_index])
      alt_index = alt_index % 2 + 1
    end
  end

  x = 0.5
  m_0 = 2 * sqrt(x)
  #h_1 = 0.4

  gamma_list2 = [1.5,2.0,2.5]   # <<< two gamma value

  for gamma in gamma_list2

    println("\nRunning for gamma = $gamma")

    h = -im * gamma
    psi0 = MPS(s, state)   # reset state for each gamma

    gates1 = ITensor[]
    gates3 = ITensor[]
    gates4 = ITensor[]

    ################# Hopping term #################
    for j in 1:2:N-2
      s1 = s[j]
      s2 = s[j+2]
      s3 = s[j+1]

      hj = (x+h) * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

      Gj = exp(-im * tau * hj)
      push!(gates1, Gj)
    end

    ################# Electric field #################
    for j in 2:2:N-1
      hj = f * op("Sz", s[j])
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates3, Gj)
    end

    ################# Mass + non-Hermitian #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)
    end

    ################# Full Trotter step #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# Output #################
    filename = "SvN_HIfullgamma$(gamma)_64.dat"
    open(filename, "w") do io
      for t in 0:tau:ttotal
        S_vN = ent1(psi0)
        write(io, "$t\t$S_vN\n")
        println("gamma=$gamma  t=$t  SvN=$S_vN")

        psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
        normalize!(psi0)
      end
    end
  end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
  ############################
  # PARAMETERS
  ############################
  f = 1
  N = 127
  tau = 0.1
  ttotal = 100
  cutoff = 1E-8
  bond_dim_max = 1000

  s = siteinds("S=1/2", N)

  # Initial staggered state
  state = String[]
  alt = ["Dn", "Up"]
  alt_index = 1
  for i in 1:N
    if iseven(i)
      push!(state, "Dn")
    else
      push!(state, alt[alt_index])
      alt_index = alt_index % 2 + 1
    end
  end

  x = 0.5
  m_0 = 2 * sqrt(x)
  #h_1 = 0.4

  gamma_list3 = [3.0,3.5,4.0]   # <<< two gamma value

  for gamma in gamma_list3

    println("\nRunning for gamma = $gamma")

    h = -im * gamma
    psi0 = MPS(s, state)   # reset state for each gamma

    gates1 = ITensor[]
    gates3 = ITensor[]
    gates4 = ITensor[]

    ################# Hopping term #################
    for j in 1:2:N-2
      s1 = s[j]
      s2 = s[j+2]
      s3 = s[j+1]

      hj = (x+h) * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

      Gj = exp(-im * tau * hj)
      push!(gates1, Gj)
    end

    ################# Electric field #################
    for j in 2:2:N-1
      hj = f * op("Sz", s[j])
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates3, Gj)
    end

    ################# Mass + non-Hermitian #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)
    end

    ################# Full Trotter step #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# Output #################
    filename = "SvN_HIfullgamma$(gamma)_64.dat"
    open(filename, "w") do io
      for t in 0:tau:ttotal
        S_vN = ent1(psi0)
        write(io, "$t\t$S_vN\n")
        println("gamma=$gamma  t=$t  SvN=$S_vN")

        psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
        normalize!(psi0)
      end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    x = 0.5
    m_0 = 2 * sqrt(x)

    gamma_list = [0.1,0.5,1]   # <<< THREE GAMMAS

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER GAMMA
    ############################
    for gamma in gamma_list

     println("\nRunning for gamma = $gamma")
     h = -im * gamma

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        if (j >= 1) && (j <= b-2)
            coeff = x + h
        else
            coeff = x
        end

        hj = coeff * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HIsub_gamma_$(gamma).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("gamma=$gamma  t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    x = 0.5
    m_0 = 2 * sqrt(x)

    gamma_list1 = [1.5,2.0,2.5]   # <<< THREE GAMMAS

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER GAMMA
    ############################
    for gamma in gamma_list1

     println("\nRunning for gamma = $gamma")
     h = -im * gamma

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        if (j >= 1) && (j <= b-2)
            coeff = x + h
        else
            coeff = x
        end

        hj = coeff * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HIsub_gamma_$(gamma).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("gamma=$gamma  t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    x = 0.5
    m_0 = 2 * sqrt(x)

    gamma_list2 = [3.0,3.5,4.0]   # <<< THREE GAMMAS

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER GAMMA
    ############################
    for gamma in gamma_list2

     println("\nRunning for gamma = $gamma")
     h = -im * gamma

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        if (j >= 1) && (j <= b-2)
            coeff = x + h
        else
            coeff = x
        end

        hj = coeff * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HIsub_gamma_$(gamma).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("gamma=$gamma  t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

N=7
b=(N ÷ 2)
for j in 1:2:N-2
 if (1 <= j <= b-2) || (b< j <= N-2)
  println(j)
 else
  println("link=$j")
 end
end

N=7
b=(N ÷ 2)
for j in 1:2:N-2
 if j <= b-2 || j >= b+1
    println(j)
 else
    println("link=$j")
 end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 80
    cutoff = 1e-8
    maxdim = 1000

    x = 0.5
    m_0 = 2 * sqrt(x)

    gamma_list1 = [0.1,0.5,1.0]   # <<< THREE GAMMAS

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER GAMMA
    ############################
    for gamma in gamma_list1

     println("\nRunning for gamma = $gamma")
     h = -im * gamma

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]
        if (1<= j <= b-2) || (b<= j <= N-2)
          hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
        else
          hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
        end
        Gj = exp(-im * (tau) *hj)
        push!(gates1, Gj)
     end

        

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HInolink_gamma_$(gamma).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("gamma=$gamma  t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 80
    cutoff = 1e-8
    maxdim = 1000

    x = 0.5
    m_0 = 2 * sqrt(x)

    gamma_list2 = [1.5,2.0,2.5]   # <<< THREE GAMMAS

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER GAMMA
    ############################
    for gamma in gamma_list2

     println("\nRunning for gamma = $gamma")
     h = -im * gamma

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]
        if (1<= j <= b-2) || (b<= j <= N-2)
          hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
        else
          hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
        end
        Gj = exp(-im * (tau) *hj)
        push!(gates1, Gj)
     end

        

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HInolink_gamma_$(gamma).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("gamma=$gamma  t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 80
    cutoff = 1e-8
    maxdim = 1000

    x = 0.5
    m_0 = 2 * sqrt(x)

    gamma_list3 = [3.0,3.5,4.0]   # <<< THREE GAMMAS

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER GAMMA
    ############################
    for gamma in gamma_list3

     println("\nRunning for gamma = $gamma")
     h = -im * gamma

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]
        if (1<= j <= b-2) || (b<= j <= N-2)
          hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
        else
          hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
        end
        Gj = exp(-im * (tau) *hj)
        push!(gates1, Gj)
     end

        

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HInolink_gamma_$(gamma).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("gamma=$gamma  t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 80
    cutoff = 1e-8
    maxdim = 1000

    x = 1.5
    m_0 = 2 * sqrt(x)

    gamma_list1 = [0.2,0.5,1.0,1.5]   # <<< THREE GAMMAS

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER GAMMA
    ############################
    for gamma in gamma_list1

     println("\nRunning for gamma = $gamma")
     h = -im * gamma

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]
        if (1<= j <= b-2) || (b<= j <= N-2)
          hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
        else
          hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
        end
        Gj = exp(-im * (tau) *hj)
        push!(gates1, Gj)
     end

        

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HInolink_gamma2_$(gamma).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("gamma=$gamma  t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 80
    cutoff = 1e-8
    maxdim = 1000

    x = 1.5
    m_0 = 2 * sqrt(x)

    gamma_list2 = [2.0,2.5,3.0,3.5]   # <<< THREE GAMMAS

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER GAMMA
    ############################
    for gamma in gamma_list2

     println("\nRunning for gamma = $gamma")
     h = -im * gamma

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]
        if (1<= j <= b-2) || (b<= j <= N-2)
          hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
        else
          hj = x*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
        end
        Gj = exp(-im * (tau) *hj)
        push!(gates1, Gj)
     end

        

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HInolink_gamma2_$(gamma).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("gamma=$gamma  t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 0.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list1 = [0.2,0.4,0.6,0.8]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list1
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        if (j >= 1) && (j <= b-2)
            coeff = x + h
        else
            coeff = x
        end

        hj = coeff * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HIsub_x_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 0.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list2 = [1.0,4.0,8.0,16.0]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list2
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        if (j >= 1) && (j <= b-2)
            coeff = x + h
        else
            coeff = x
        end

        hj = coeff * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HIsub_x_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 0.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list1 = [0.2,0.4,0.6,0.8]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list1
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = (x+h)*(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HIfull_x_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 0.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list2 = [1,4,8,16]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list2
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = (x+h)*(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HIfull_x_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 127
  b=(N/2)
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.5
  h= -im*h_1
  
  x    = 2.0
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x+h)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HIfull_x_2.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 1.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list1 = [0.2,0.4,0.6,0.8]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list1
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = (x+h)*(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HIfull2_x_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 1.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list2 = [1,2,4,8]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list2
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = (x+h)*(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HIfull2_x_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 40
    cutoff = 1e-8
    maxdim = 1000

    gamma = 1.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)
    x = 2 
    m_0 = 2 * sqrt(x)
    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    #for x in x_list2
     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = (x+h)*(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    #filename = "SvN_HIfull2_x_$(x).dat"
    #open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            #write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    #end
  #end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 40
    cutoff = 1e-8
    maxdim = 1000

    gamma = 1.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)
    x = 4 
    m_0 = 2 * sqrt(x)
    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    #for x in x_list2
     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = (x+h)*(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = f * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    #filename = "SvN_HIfull2_x_$(x).dat"
    #open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            #write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    #end
  #end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
  ############################
  # PARAMETERS
  ############################
  f = 1
  N = 127
  tau = 0.1
  ttotal = 100
  cutoff = 1E-8
  bond_dim_max = 1000

  s = siteinds("S=1/2", N)

  # Initial staggered state
  state = String[]
  alt = ["Dn", "Up"]
  alt_index = 1
  for i in 1:N
    if iseven(i)
      push!(state, "Dn")
    else
      push!(state, alt[alt_index])
      alt_index = alt_index % 2 + 1
    end
  end

  x = 0.5
  m_0 = 2 * sqrt(x)
  #h_1 = 0.4

  gamma_list1 = [0.1,0.5,1]   # <<< two gamma value

  for gamma in gamma_list1

    println("\nRunning for gamma = $gamma")

    h = -im * gamma
    psi0 = MPS(s, state)   # reset state for each gamma

    gates1 = ITensor[]
    gates3 = ITensor[]
    gates4 = ITensor[]

    ################# Hopping term #################
    for j in 1:2:N-2
      s1 = s[j]
      s2 = s[j+2]
      s3 = s[j+1]

      hj = x * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

      Gj = exp(-im * tau * hj)
      push!(gates1, Gj)
    end

    ################# Electric field #################
    for j in 2:2:N-1
      hj = (f) * op("Sz", s[j])
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates3, Gj)
    end

    ################# Mass + non-Hermitian #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0+h) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0+h) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)
    end

    ################# Full Trotter step #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# Output #################
    filename = "SvN_numopgamma$(gamma)_64.dat"
    open(filename, "w") do io
      for t in 0:tau:ttotal
        S_vN = ent1(psi0)
        write(io, "$t\t$S_vN\n")
        println("gamma=$gamma  t=$t  SvN=$S_vN")

        psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
        normalize!(psi0)
      end
    end
  end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
  ############################
  # PARAMETERS
  ############################
  f = 1
  N = 127
  tau = 0.1
  ttotal = 100
  cutoff = 1E-8
  bond_dim_max = 1000

  s = siteinds("S=1/2", N)

  # Initial staggered state
  state = String[]
  alt = ["Dn", "Up"]
  alt_index = 1
  for i in 1:N
    if iseven(i)
      push!(state, "Dn")
    else
      push!(state, alt[alt_index])
      alt_index = alt_index % 2 + 1
    end
  end

  x = 0.5
  m_0 = 2 * sqrt(x)
  #h_1 = 0.4

  gamma_list2 = [1.5,2,2.5]   # <<< two gamma value

  for gamma in gamma_list2

    println("\nRunning for gamma = $gamma")

    h = -im * gamma
    psi0 = MPS(s, state)   # reset state for each gamma

    gates1 = ITensor[]
    gates3 = ITensor[]
    gates4 = ITensor[]

    ################# Hopping term #################
    for j in 1:2:N-2
      s1 = s[j]
      s2 = s[j+2]
      s3 = s[j+1]

      hj = x * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

      Gj = exp(-im * tau * hj)
      push!(gates1, Gj)
    end

    ################# Electric field #################
    for j in 2:2:N-1
      hj = (f) * op("Sz", s[j])
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates3, Gj)
    end

    ################# Mass + non-Hermitian #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0+h) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0+h) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)
    end

    ################# Full Trotter step #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# Output #################
    filename = "SvN_numopgamma$(gamma)_64.dat"
    open(filename, "w") do io
      for t in 0:tau:ttotal
        S_vN = ent1(psi0)
        write(io, "$t\t$S_vN\n")
        println("gamma=$gamma  t=$t  SvN=$S_vN")

        psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
        normalize!(psi0)
      end
    end
  end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
  ############################
  # PARAMETERS
  ############################
  f = 1
  N = 127
  tau = 0.1
  ttotal = 100
  cutoff = 1E-8
  bond_dim_max = 1000

  s = siteinds("S=1/2", N)

  # Initial staggered state
  state = String[]
  alt = ["Dn", "Up"]
  alt_index = 1
  for i in 1:N
    if iseven(i)
      push!(state, "Dn")
    else
      push!(state, alt[alt_index])
      alt_index = alt_index % 2 + 1
    end
  end

  x = 0.5
  m_0 = 2 * sqrt(x)
  #h_1 = 0.4

  gamma_list3 = [3,3.5,4]   # <<< two gamma value

  for gamma in gamma_list3

    println("\nRunning for gamma = $gamma")

    h = -im * gamma
    psi0 = MPS(s, state)   # reset state for each gamma

    gates1 = ITensor[]
    gates3 = ITensor[]
    gates4 = ITensor[]

    ################# Hopping term #################
    for j in 1:2:N-2
      s1 = s[j]
      s2 = s[j+2]
      s3 = s[j+1]

      hj = x * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

      Gj = exp(-im * tau * hj)
      push!(gates1, Gj)
    end

    ################# Electric field #################
    for j in 2:2:N-1
      hj = (f) * op("Sz", s[j])
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates3, Gj)
    end

    ################# Mass + non-Hermitian #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0+h) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0+h) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)
    end

    ################# Full Trotter step #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# Output #################
    filename = "SvN_numopgamma$(gamma)_64.dat"
    open(filename, "w") do io
      for t in 0:tau:ttotal
        S_vN = ent1(psi0)
        write(io, "$t\t$S_vN\n")
        println("gamma=$gamma  t=$t  SvN=$S_vN")

        psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
        normalize!(psi0)
      end
    end
  end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
  ############################
  # PARAMETERS
  ############################
  f = 1
  N = 127
  tau = 0.1
  ttotal = 100
  cutoff = 1E-8
  bond_dim_max = 1000

  s = siteinds("S=1/2", N)

  # Initial staggered state
  state = String[]
  alt = ["Dn", "Up"]
  alt_index = 1
  for i in 1:N
    if iseven(i)
      push!(state, "Dn")
    else
      push!(state, alt[alt_index])
      alt_index = alt_index % 2 + 1
    end
  end

  x = 0.5
  m_0 = 2 * sqrt(x)
  #h_1 = 0.4

  gamma_list_1 = [0.1,0.5,1]   # <<< two gamma value

  for gamma in gamma_list_1

    println("\nRunning for gamma = $gamma")

    h = -im * gamma
    psi0 = MPS(s, state)   # reset state for each gamma

    gates1 = ITensor[]
    gates3 = ITensor[]
    gates4 = ITensor[]

    ################# Hopping term #################
    for j in 1:2:N-2
      s1 = s[j]
      s2 = s[j+2]
      s3 = s[j+1]

      hj = x * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

      Gj = exp(-im * tau * hj)
      push!(gates1, Gj)
    end

    ################# Electric field #################
    for j in 2:2:N-1
      hj = (f+h) * op("Sz", s[j])
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates3, Gj)
    end

    ################# Mass + non-Hermitian #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)
    end

    ################# Full Trotter step #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# Output #################
    filename = "SvN_HEgamma$(gamma)_64.dat"
    open(filename, "w") do io
      for t in 0:tau:ttotal
        S_vN = ent1(psi0)
        write(io, "$t\t$S_vN\n")
        println("gamma=$gamma  t=$t  SvN=$S_vN")

        psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
        normalize!(psi0)
      end
    end
  end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
  ############################
  # PARAMETERS
  ############################
  f = 1
  N = 127
  tau = 0.1
  ttotal = 100
  cutoff = 1E-8
  bond_dim_max = 1000

  s = siteinds("S=1/2", N)

  # Initial staggered state
  state = String[]
  alt = ["Dn", "Up"]
  alt_index = 1
  for i in 1:N
    if iseven(i)
      push!(state, "Dn")
    else
      push!(state, alt[alt_index])
      alt_index = alt_index % 2 + 1
    end
  end

  x = 0.5
  m_0 = 2 * sqrt(x)
  #h_1 = 0.4

  gamma_list_2 = [1.5,2.0,2.5]   # <<< three gamma value

  for gamma in gamma_list_2

    println("\nRunning for gamma = $gamma")

    h = -im * gamma
    psi0 = MPS(s, state)   # reset state for each gamma

    gates1 = ITensor[]
    gates3 = ITensor[]
    gates4 = ITensor[]

    ################# Hopping term #################
    for j in 1:2:N-2
      s1 = s[j]
      s2 = s[j+2]
      s3 = s[j+1]

      hj = x * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

      Gj = exp(-im * tau * hj)
      push!(gates1, Gj)
    end

    ################# Electric field #################
    for j in 2:2:N-1
      hj = (f+h) * op("Sz", s[j])
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates3, Gj)
    end

    ################# Mass + non-Hermitian #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)
    end

    ################# Full Trotter step #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# Output #################
    filename = "SvN_HEgamma$(gamma)_64.dat"
    open(filename, "w") do io
      for t in 0:tau:ttotal
        S_vN = ent1(psi0)
        write(io, "$t\t$S_vN\n")
        println("gamma=$gamma  t=$t  SvN=$S_vN")

        psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
        normalize!(psi0)
      end
    end
  end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
  ############################
  # PARAMETERS
  ############################
  f = 1
  N = 127
  tau = 0.1
  ttotal = 100
  cutoff = 1E-8
  bond_dim_max = 1000

  s = siteinds("S=1/2", N)

  # Initial staggered state
  state = String[]
  alt = ["Dn", "Up"]
  alt_index = 1
  for i in 1:N
    if iseven(i)
      push!(state, "Dn")
    else
      push!(state, alt[alt_index])
      alt_index = alt_index % 2 + 1
    end
  end

  x = 0.5
  m_0 = 2 * sqrt(x)
  #h_1 = 0.4

  gamma_list_3 = [3.0,3.5,4.0]   # <<< three gamma value

  for gamma in gamma_list_3

    println("\nRunning for gamma = $gamma")

    h = -im * gamma
    psi0 = MPS(s, state)   # reset state for each gamma

    gates1 = ITensor[]
    gates3 = ITensor[]
    gates4 = ITensor[]

    ################# Hopping term #################
    for j in 1:2:N-2
      s1 = s[j]
      s2 = s[j+2]
      s3 = s[j+1]

      hj = x * (op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

      Gj = exp(-im * tau * hj)
      push!(gates1, Gj)
    end

    ################# Electric field #################
    for j in 2:2:N-1
      hj = (f+h) * op("Sz", s[j])
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates3, Gj)
    end

    ################# Mass + non-Hermitian #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)
    end

    ################# Full Trotter step #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# Output #################
    filename = "SvN_HEgamma$(gamma)_64.dat"
    open(filename, "w") do io
      for t in 0:tau:ttotal
        S_vN = ent1(psi0)
        write(io, "$t\t$S_vN\n")
        println("gamma=$gamma  t=$t  SvN=$S_vN")

        psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
        normalize!(psi0)
      end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 0.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list1 = [0.2,0.4,0.6,0.8]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list1
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = x *(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = (f+h) * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HE_x_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 0.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list2 = [1,1.2,1.4,1.6]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list2
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = x *(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = (f+h) * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HE_x_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 0.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list3 = [1.8,2.0,2.2,2.4]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list3
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = x *(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = (f+h) * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HE_x_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 1.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list1 = [0.2,0.4,0.6,0.8]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list1
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = x *(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = (f+h) * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HE_x2_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 1.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list1 = [0.2,0.4,0.6,0.8]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list1
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = x *(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = (f+h) * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,(gates5))

    ################# OUTPUT #################
    filename = "SvN_HE_x2NR_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 1.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list2 = [1,1.2,1.4,1.6]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list2
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = x *(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = (f+h) * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HE_x2_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 1.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list2 = [1,1.2,1.4,1.6]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list2
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = x *(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = (f+h) * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,(gates5))

    ################# OUTPUT #################
    filename = "SvN_HE_x2NR_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N ÷ 2
  orthogonalize!(psi, b)
  U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
  for n = 1:dim(S, 1)
    p = S[n, n]^2
    SvN -= p * log(p)
  end
  return SvN
end

let
    ############################
    # PARAMETERS
    ############################
    f = 1.0
    N = 127
    b = N ÷ 2
    tau = 0.1
    ttotal = 100
    cutoff = 1e-8
    maxdim = 1000

    gamma = 1.5
    h = -im * gamma
    #m_0 = 2 * sqrt(x)

    x_list3 = [1.8,2.0,2.2,2.4]   # <<< 4 x-value

    ############################
    # SITE INDICES
    ############################
    s = siteinds("S=1/2", N)

    ############################
    # INITIAL STATE
    ############################
    state = String[]
    alt = ["Dn", "Up"]
    alt_index = 1
    for i in 1:N
        if iseven(i)
            push!(state, "Dn")
        else
            push!(state, alt[alt_index])
            alt_index = alt_index % 2 + 1
        end
    end

    ############################
    # LOOP OVER X
    ############################
    for x in x_list3
     m_0 = 2 * sqrt(x)

     println("\nRunning for x = $x,m_0= $m_0")

     psi0 = MPS(s, state)

     gates1 = ITensor[]
     gates3 = ITensor[]
     gates4 = ITensor[]

     ################# HOPPING #################
     for j in 1:2:N-2
        s1, s2, s3 = s[j], s[j+2], s[j+1]

        hj = x *(op("S+", s1) * op("Sx", s3) * op("S-", s2) +op("S+", s2) * op("Sx", s3) * op("S-", s1))

        push!(gates1, exp(-im * tau * hj))
        end

     ################# ELECTRIC FIELD #################
     for j in 2:2:N-1
        hj = (f+h) * op("Sz", s[j])
        push!(gates3, exp(-im * tau/2 * hj))
     end

    ################# MASS #################
    p = 0
    for j in 1:2:N
      if p == 0
        hj = (m_0) * op("Sz", s[j])
        p = 1
      else
        hj = -(m_0) * op("Sz", s[j])
        p = 0
      end
      Gj = exp(-im * (tau / 2) * hj)
      push!(gates4, Gj)  
    end

    ################# FULL TROTTER STEP #################
    gates5=append!(gates3,gates4)
    gates=append!(gates5,gates1,reverse(gates5))

    ################# OUTPUT #################
    filename = "SvN_HE_x2_$(x).dat"
    open(filename, "w") do io
        for t in 0:tau:ttotal
            S_vN = ent1(psi0)
            write(io, "$t\t$S_vN\n")
            println("x=$x t=$t  SvN=$S_vN")

            psi0 = apply(gates, psi0; cutoff, maxdim)
            normalize!(psi0)
        end
    end
  end
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 11
  b=(N/2)
  tau = 0.08
  ttotal = 10
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 0
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = (m_0+h)*op("Sz",s1)
     p=1
    elseif p==1
     hj = -(m_0+h)*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_numopx0.5_6.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 11
  b=(N/2)
  tau = 0.08
  ttotal = 10
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 0
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau/2) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = (m_0+h)*op("Sz",s1)
     p=1
    elseif p==1
     hj = -(m_0+h)*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates1,gates5,reverse(gates1))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_numopx0.5interchange_6.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 11
  b=(N/2)
  tau = 0.08
  ttotal = 10
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 0
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = (m_0+h)*op("Sz",s1)
     p=1
    elseif p==1
     hj = -(m_0+h)*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_numopx0.5NR_6.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 11
  b=(N/2)
  tau = 0.08
  ttotal = 10
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 0
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau/2) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = f*op("Sz",s1)
    Gj = exp(-im * (tau) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = (m_0+h)*op("Sz",s1)
     p=1
    elseif p==1
     hj = -(m_0+h)*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates1,gates5,(gates1))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_numopx0.5NRinterchange_6.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff1 = 1E-10

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HEgammacutoffcheck_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff=cutoff1,maxdim = bond_dim_max)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N=127
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.4
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]
        
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
##################Gauss#########
################################

################################  

#  for j in 1:2:N-1
#   for j in 1:N-1
#    s1 = s[j]
#    s2 = s[j+1]
#    hj = 0.5*gamma*op("Sx",s1)*op("Sx",s2)

  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  
   
  
  
  
  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max1=2000
  filename = "SvN_HEgammabonddimcheck_64.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
    S_vN= ent1(psi0)
    psi0 = apply(gates,psi0;cutoff=cutoff,maxdim = bond_dim_max1)
    normalize!(psi0)
    write(io, "$t \t $S_vN \n")
    println("$t $S_vN")   
   end
  end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 127
  b=(N/2)
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1= 0.5
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HE_time.dat"
  open(filename, "w") do io
   for t in 0:tau:ttotal
     t_start = time_ns()
     S_vN = ent1(psi0)
     t_end  = time_ns()
     EE_time = (t_end - t_start) * 1e-9  # seconds

     write(io, "$t\t$S_vN\t$EE_time\n")

     println("gamma=$h  t=$t  SvN=$S_vN  EE_time=$(round(EE_time, digits=6)) s")

     psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
     normalize!(psi0)
  end
 end   
return
end

using ITensors, ITensorMPS

mutable struct SizeObserver <: AbstractObserver
end

function ent1(psi::MPS)
  N = length(psi)
  SvN = 0.0
  b = N÷2
  #b=65
  orthogonalize!(psi,b)
  U,S,V = svd(psi[b],(linkind(psi,b-1),siteind(psi,b)))
  D = dim(S,1)
  for n=1:D
    p = S[n,n]^2
    SvN -= p*log(p)
  end
  return SvN
end

using ITensorMPS, ITensors

mutable struct SizeObserver <: AbstractObserver
end

function ITensorMPS.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end


let
  f=1
  N= 127
  b=(N/2)
  tau = 0.1
  ttotal = 100
  s = siteinds("S=1/2", N)
  # Make gates (1,2),(2,3),(3,4),...
  gates1  = ITensor[]
  gates2  = ITensor[]
  gates3  = ITensor[]
  gates4  = ITensor[]
  gates5  = ITensor[]
  gates6  = ITensor[]
  gates   = ITensor[]
  p=0
  
  state=[]
  
  alt = ["Dn", "Up"]
  alt_index = 1

  for i in 1:N
    if iseven(i)
        push!(state, "Dn")
    else
        push!(state, alt[alt_index])
        alt_index = alt_index % 2 + 1
    end
  end
  psi0 = MPS(s, state)
  h_1=0.5
  h= -im*h_1
  
  x    = 0.5
  m_0  = 2*(x)^(0.5)
  @show x
  @show m_0
  
  #gamma= 0.75
  #g1=0.4
  #J=0.4
  
  
#################hermitian term+hopping term######################
   
  for j in 1:2:N-2
    s1 = s[j]
    s2 = s[j+2]
    s3 = s[j+1]   
    hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    Gj = exp(-im * (tau) *hj)
    push!(gates1, Gj)
  end

  #for j in 3:4:N-2
    #s1 = s[j]
    #s2 = s[j+2]
    #s3 = s[j+1]
   
    #hj = (n_0*op("Sx", s1)*op("Sx", s2) + n_0*op("Sy", s1)*op("Sy", s2))*op("Sz",s3)
    #hj = (x)*(op("S+",s1)*op("Sx",s3)*op("S-",s2)+op("S+",s2)*op("Sx",s3)*op("S-",s1))
    #Gj = exp(-im * (tau) *hj)
    #push!(gates2, Gj)
  #end  

    
 ############################################### 
  
  
  for j in 2:2:N-1
    s1 = s[j]
    hj = (f+h)*op("Sz",s1)
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates3, Gj)
  end  


######################mass###################  
  p=0 
  for j in 1:2:N
    s1 = s[j]
    if p==0
     hj = m_0*op("Sz",s1)
     p=1
    elseif p==1
     hj = -m_0*op("Sz",s1) 
     p=0
    end 
    Gj = exp(-im * (tau / 2) *hj)
    push!(gates4, Gj)
  end
  

  gates5=append!(gates3,gates4)
  gates=append!(gates5,gates1,reverse(gates5))

  #@show psi0
  cutoff = 1E-8

  obs = SizeObserver()
  #y_list=[] #x_list=[]
  bond_dim_max=1000
  filename = "SvN_HE_time_chi.dat"

  open(filename, "w") do io
   write(io, "# t\tSvN\tEE_time(s)\tchi_mid\tchi_max\n")

  for t in 0:tau:ttotal

    # ---------- EE timing ----------
    t_start = time_ns()
    S_vN = ent1(psi0)
    EE_time = (time_ns() - t_start) * 1e-9

    # ---------- bond dimensions ----------
    b = length(psi0) ÷ 2
    chi_mid = dim(linkind(psi0, b))
    chi_max = maxlinkdim(psi0)

    # ---------- output ----------
    write(io, "$t\t$S_vN\t$EE_time\t$chi_mid\t$chi_max\n")

    println(
      "t=$t  SvN=$(round(S_vN,digits=6))  " *
      "EE_time=$(round(EE_time,digits=6)) s  " *
      "χ_mid=$chi_mid  χ_max=$chi_max"
    )

    # ---------- time evolution ----------
    psi0 = apply(gates, psi0; cutoff, maxdim=bond_dim_max)
    normalize!(psi0)
  end
 end
return
end


