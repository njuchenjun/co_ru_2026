subroutine pes_ru12co(x,e,dedx)
  use nnmod_pes1, only: init_pes1,evaluate_pes1,deallocate_pes1
  implicit none
  real*8,intent(in) :: x(3,14)
  real*8,intent(out) :: e,dedx(3,14)
  real*8 :: c(3,22),v,dvdc(3,22)
  character(kind=1,len=2),parameter :: elems(22)=["Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","C ","O "]
  integer,save :: ifirst = 0
  if(ifirst .eq. 0) then
    call init_pes1
    ifirst = 1
  end if
  
  c(1:3, 1)=[   0.79107d0,   2.43599d0,   0.00000d0 ]
  c(1:3, 2)=[   3.16424d0,   1.06584d0,   0.00000d0 ]
  c(1:3, 3)=[   0.79107d0,  -0.30431d0,   0.00000d0 ]
  c(1:3, 4)=[   3.16424d0,  -1.67446d0,   0.00000d0 ]
  c(1:3, 5)=[   1.58209d0,   1.06584d0,   2.13254d0 ]
  c(1:3, 6)=[   3.95526d0,  -0.30431d0,   2.13254d0 ]
  c(1:3, 7)=[   1.58209d0,   3.80613d0,   2.13254d0 ]
  c(1:3, 8)=[   3.95526d0,   2.43599d0,   2.13254d0 ]
 !c(1:3, 9)=[   0.79107d0,   2.43599d0,   4.26507d0 ]
 !c(1:3,10)=[   3.16424d0,   1.06584d0,   4.26507d0 ]
 !c(1:3,11)=[   0.79107d0,  -0.30431d0,   4.26507d0 ]
 !c(1:3,12)=[   3.16424d0,  -1.67446d0,   4.26507d0 ]
 !c(1:3,13)=[   1.58176d0,   1.06581d0,   6.41455d0 ]
 !c(1:3,14)=[   3.95493d0,  -0.30434d0,   6.41455d0 ]
 !c(1:3,15)=[   1.58176d0,   3.80611d0,   6.41455d0 ]
 !c(1:3,16)=[   3.95493d0,   2.43596d0,   6.41455d0 ]
 !c(1:3,17)=[   0.79117d0,   2.43599d0,   8.47202d0 ]
 !c(1:3,18)=[   3.16433d0,   1.06584d0,   8.47202d0 ]
 !c(1:3,19)=[   0.79117d0,  -0.30431d0,   8.47202d0 ]
 !c(1:3,10)=[   3.16433d0,  -1.67446d0,   8.47202d0 ]
 !c(1:3,21)=[   1.15471d0,   0.28944d0,  16.29591d0 ]
 !c(1:3,22)=[   1.05598d0,  -0.05304d0,  15.26985d0 ]
  c(1:3,9:22)=x
  call evaluate_pes1(0,1,22,elems,c,v,dvdc)
  e = v
  dedx = dvdc(1:3,9:22)
  return
end subroutine

subroutine pes_ru8co(x,e,dedx)
  use nnmod_pes1, only: init_pes1,evaluate_pes1,deallocate_pes1
  implicit none
  real*8,intent(in) :: x(3,10)
  real*8,intent(out) :: e,dedx(3,10)
  real*8 :: c(3,22),v,dvdc(3,22)
  character(kind = 1,len = 2),parameter :: elems(22)=["Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","Ru","C ","O "]
  integer,save :: ifirst = 0
  if (ifirst .eq. 0) then
     call init_pes1
     ifirst = 1
  end if
  c(1:3, 1)=[   0.79107d0,   2.43599d0,   0.00000d0 ]
  c(1:3, 2)=[   3.16424d0,   1.06584d0,   0.00000d0 ]
  c(1:3, 3)=[   0.79107d0,  -0.30431d0,   0.00000d0 ]
  c(1:3, 4)=[   3.16424d0,  -1.67446d0,   0.00000d0 ]
  c(1:3, 5)=[   1.58209d0,   1.06584d0,   2.13254d0 ]
  c(1:3, 6)=[   3.95526d0,  -0.30431d0,   2.13254d0 ]
  c(1:3, 7)=[   1.58209d0,   3.80613d0,   2.13254d0 ]
  c(1:3, 8)=[   3.95526d0,   2.43599d0,   2.13254d0 ]
  c(1:3, 9)=[   0.79107d0,   2.43599d0,   4.26507d0 ]
  c(1:3,10)=[   3.16424d0,   1.06584d0,   4.26507d0 ]
  c(1:3,11)=[   0.79107d0,  -0.30431d0,   4.26507d0 ]
  c(1:3,12)=[   3.16424d0,  -1.67446d0,   4.26507d0 ]
 !c(1:3,13)=[   1.58176d0,   1.06581d0,   6.41455d0 ]
 !c(1:3,14)=[   3.95493d0,  -0.30434d0,   6.41455d0 ]
 !c(1:3,15)=[   1.58176d0,   3.80611d0,   6.41455d0 ]
 !c(1:3,16)=[   3.95493d0,   2.43596d0,   6.41455d0 ]
 !c(1:3,17)=[   0.79117d0,   2.43599d0,   8.47202d0 ]
 !c(1:3,18)=[   3.16433d0,   1.06584d0,   8.47202d0 ]
 !c(1:3,19)=[   0.79117d0,  -0.30431d0,   8.47202d0 ]
 !c(1:3,20)=[   3.16433d0,  -1.67446d0,   8.47202d0 ]
 !c(1:3,21)=[   1.15471d0,   0.28944d0,  16.29591d0 ]
 !c(1:3,22)=[   1.05598d0,  -0.05304d0,  15.26985d0 ]
  c(1:3,13:22)=x
  call evaluate_pes1(0,1,22,elems,c,v,dvdc)
  e = v
  dedx = dvdc(1:3,13:22)
  return
end subroutine

subroutine pes_ruFixco(x,e,dedx)
  implicit none
  real*8,intent(in) :: x(3,2)
  real*8,intent(out) :: e,dedx(3,2)
  call co_ru_pes(x,e,1,dedx)
  dedx = -dedx
  return
end subroutine
