module nnmod_pes1
     integer(kind=4),parameter :: intype=4,retype=8,atomdim=3
     real(kind=retype) :: pi=dacos(-1.d0)
     character(kind=1,len=80) :: pesdir='./pes8.dat'
     integer(kind=intype) :: pesid=1
     include "pes_hd.inc1.f90"
end module

