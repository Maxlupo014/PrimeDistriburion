module facto
contains
  function fact(a) result(n)
    implicit none
    integer, parameter :: ik= selected_int_kind(15)
    integer(kind=ik) :: a, i
    real::n
    n=1
    if(a==0)then
      n=1
    else
      do i=1, a
        n=i*n
      end do
    end if
    end function fact
end module facto


program primi
use facto
implicit none
integer, parameter :: rk1= selected_int_kind(6)
integer, parameter :: rk2= selected_int_kind(15)
integer, parameter :: ik= selected_int_kind(15)
integer, parameter :: rk= selected_real_kind(15)
integer(kind=ik), dimension(:), allocatable :: primo
real(kind=rk) :: med, dev=0, med1=0, dev1=0, mid, medbin, X2, X2p, stat
real(kind=rk2), dimension(:), allocatable :: d, delc, ext, c
integer(kind=ik) :: i,r,m, w=1, l, j=0, k, n=40, max=1000000, rag=1000, counter, pokt=4, n0, n1
integer(kind=ik), dimension(16) :: arr
integer(kind=ik), dimension(:), allocatable :: sums
integer(kind=ik), dimension(:), allocatable :: fin, tr, pok
real, dimension(:), allocatable :: ten, devst, logst
open(Unit=12, file='primes200.txt')


allocate(primo(max))
allocate(d(max-2))
allocate(c(max-2))

do i=1, max
  read(12,*) primo(i)
end do

allocate(ext((primo(max)-primo(1))/n+1))
allocate(sums((primo(max)-primo(1))/n+1))
allocate(fin(n+1))
fin=0
sums=0
arr=0

do i=1, max
!𝑛(ln(𝑛)+ln(ln(𝑛)−1)
l=((primo(i)-primo(1))/n)
!l=(((primo(i)*(log(real(primo(i))+log(log(real(primo(i)-1)))))))-((primo(1)*(log(real(primo(1))+log(log(real(primo(1)-1))))))))/n
    if (l==w) then
        sums(w)=sums(w)+1
    else
        sums(l)=1
    end if
w=l
end do

do i=1, (primo(max)-primo(1))/n+1
  k=sums(i)+1
  fin(k)=fin(k)+1
end do

med=real(sum(sums))/((primo(max)-primo(1))/n+1)
print*, med

do i=1, (primo(max)-primo(1))/n+1

    dev=dev+(real(sums(i))-med)**2
  !  print*, dev

end do

dev=dev/((primo(max)-primo(1))/n)

!--------------logaritmo

do i=1, (primo(max)-primo(1))/n+1
  ext(i)= real(n)/(log(real(primo(1)+(n+1)*i)))
  !print*, sums(i), ext(i)
end do


!-------------------------calcolo delta
do i=2, max-1

d(i)=(real(primo(i+1)-primo(i)))/(real(primo(i)-primo(i-1)))
c(i)=(real(primo(i+1)-primo(i)))
end do

do i=1, max-2
write(13,*)i, c(i)
end do

do i=2, max-1
    if((d(i)<1))then
      med1=med1-1./d(i)
    else
      med1=med1+d(i)-1
    end if
end do

med1=med1/(max-2)
dev1=0
do i=2, max-1

    if((d(i)<1))then
      dev1=dev1+(-1./d(i)-med1)**2
    else
      dev1=dev1+(d(i)-med1-1)**2
    end if

end do

dev1=dev1/(max-3)

do i=2, max-1
  if((d(i)<8).and.(-1./d(i)>-8))then
    if((d(i)<1))then
      m=-1./d(i)
    else
      m=d(i)-1
    end if
    !print*, m
    do r=1, 16
    !  print*, r
      if((m>=(r-9)) .and. (m<=(r-8)))then
      !  write(7,*) r-(8.5), arr(r)
        arr(r)=arr(r)+1
      end if

    end do

  end if
  end do
  !print*, dev1, med1
  !print*, med1 ma che cass
  do r=1, 16

    write(7,*) r-(8.5), arr(r), ((max-2)*16./15)*(1./sqrt(3.1415*2*dev1))*exp((-((real(r-8.5))-med1)**2)/(2*dev1))
  end do
    !p 'fort.7' w boxes, 'fort.7 u 1:3 w l
!-------------------------------binario
allocate(tr(max-2))
do i=1, max-2

  if(d(i)<0.99)then
    tr(i)=0
  else if(d(i)>1.01) then
    tr(i)=1
  else
    tr(i)=2
  end if
end do

medbin=0
counter=0

 do i=1, max-2
   if(tr(i)/=2)then
     medbin=medbin+tr(i)
     counter=counter+1
   end if
end do

 n0=counter-medbin
 n1=medbin

 medbin=medbin/counter
 print*, medbin, "qui!"



print*, n0, n1, counter, "helo"
allocate(pok(16))
pok=0

do i=1, max-6
if(tr(i/4)==0 .and. tr(i/4+1)==0 .and. tr(i/4+2)==0 .and. tr(i/4+3)==0 )then
  pok(1)=pok(1)+1
else if(tr(i/4)==0 .and. tr(i/4+1)==0 .and. tr(i/4+2)==0 .and. tr(i/4+3)==1 )then
  pok(2)=pok(2)+1
else if(tr(i/4)==0 .and. tr(i/4+1)==0 .and. tr(i/4+2)==1 .and. tr(i/4+3)==0 )then
  pok(3)=pok(3)+1
else if(tr(i/4)==0 .and. tr(i/4+1)==0 .and. tr(i/4+2)==1 .and. tr(i/4+3)==1 )then
  pok(4)=pok(4)+1
else if(tr(i/4)==0 .and. tr(i/4+1)==1 .and. tr(i/4+2)==0 .and. tr(i/4+3)==0 )then
  pok(5)=pok(5)+1
else if(tr(i/4)==0 .and. tr(i/4+1)==1 .and. tr(i/4+2)==0 .and. tr(i/4+3)==1 )then
  pok(6)=pok(6)+1
else if(tr(i/4)==0 .and. tr(i/4+1)==1 .and. tr(i/4+2)==1 .and. tr(i/4+3)==0 )then
  pok(7)=pok(7)+1
else if(tr(i/4)==0 .and. tr(i/4+1)==1 .and. tr(i/4+2)==1 .and. tr(i/4+3)==1 )then
  pok(8)=pok(8)+1
else if(tr(i/4)==1 .and. tr(i/4+1)==0 .and. tr(i/4+2)==0 .and. tr(i/4+3)==0 )then
  pok(9)=pok(9)+1
else if(tr(i/4)==1 .and. tr(i/4+1)==0 .and. tr(i/4+2)==0 .and. tr(i/4+3)==1 )then
  pok(10)=pok(10)+1
else if(tr(i/4)==1 .and. tr(i/4+1)==0 .and. tr(i/4+2)==1 .and. tr(i/4+3)==0 )then
  pok(11)=pok(11)+1
else if(tr(i/4)==1 .and. tr(i/4+1)==0 .and. tr(i/4+2)==1 .and. tr(i/4+3)==1 )then
  pok(12)=pok(12)+1
else if(tr(i/4)==1 .and. tr(i/4+1)==1 .and. tr(i/4+2)==0 .and. tr(i/4+3)==0 )then
  pok(13)=pok(13)+1
else if(tr(i/4)==1 .and. tr(i/4+1)==1 .and. tr(i/4+2)==0 .and. tr(i/4+3)==1 )then
  pok(14)=pok(14)+1
else if(tr(i/4)==1 .and. tr(i/4+1)==1 .and. tr(i/4+2)==1 .and. tr(i/4+3)==0 )then
  pok(15)=pok(15)+1
else if(tr(i/4)==1 .and. tr(i/4+1)==1 .and. tr(i/4+2)==1 .and. tr(i/4+3)==1 )then
  pok(16)=pok(16)+1
end if
end do

do i=1, 16
  write(21,*) i-1, pok(i)
end do

!-----------------------chi squared
X2=0
X2=(real((n1-n0)**2))/(n0+n1)!n1+n0=N => un grado di libertà
print*, X2, "chiquadrodela"
X2p=0
do i=1, 16
  X2p=X2p+(16./(max-2))*pok(i)**2
end do
X2p=X2p-max-2
print*, X2p, "chiquadropoker"

!-----------------------

!print*,med, "media"

print*, med, dev, primo(1), primo(max), "perstatist"
do i=1, n
  if(fin(i)>0.1)then
    write(4,*)i-1, real(fin(i)), ((primo(max)-primo(1))/n+1)*((med**(i-1))*exp(-med))/((1./i)*fact(i)), &
     ((primo(max)-primo(1))/n+1)*(1./sqrt(3.1415*2*dev))*exp((-((real(i-1))-med)**2)/(2*dev)), &
     3*sqrt(real(fin(i))), ((((primo(max)-primo(1))/n+1)*((med**(i-1))*exp(-med))/((1./i)*fact(i)))+ &
     (((primo(max)-primo(1))/n+1)*(1./sqrt(3.1415*2*dev))*exp((-((real(i-1))-med)**2)/(2*dev))))/2

  end if
end do

do i=1, 1000*n
  if(fin(((i/(1000))+1))>1)then
      write(9,*) real((real(i)/(1000))), &
      ((primo(max)-primo(1))/n+1)*(1./sqrt(3.1415*2*dev))*exp((-((real(real(i)/(1000)))-med)**2)/(2*dev))
  end if
end do


!------------------raggruppamento misure per fit
allocate(ten((((primo(max)-primo(1))/n+1))/rag))
allocate(devst((((primo(max)-primo(1))/n+1))/rag))
allocate(logst((((primo(max)-primo(1))/n+1))/rag))
ten=0
mid=0
do i=1, (primo(max)-primo(1))/n+1
  mid=mid+sums(i)
  if(mod(i,rag)==0)then
    ten(i/rag)=real(mid)/rag
    mid=0
  end if
end do

devst=0


do i=1, ((((primo(max)-primo(1))/n+1))/rag)
  logst(i)=real(n)/(log(real(primo(1)+(n+1)*(i+0.5)*rag)))
end do
do i=1, (primo(max)-primo(1))/n+1
  mid=mid+(sums(i)-ten(i/rag))**2
  !mid=mid+(sums(i)-logst(i/rag))**2
  if(mod(i,rag)==0)then
    devst(i/rag)=mid/(rag)
    mid=0
  end if
end do
!print*, ten
stat=0
!ogni 10 misure calcolo il valore aspettato, e vedo la media di quanto si discosta
do i=1, ((((primo(max)-primo(1))/n+1))/rag)
  stat=stat+(ten(i)-real(n)/(log(real(primo(1)+(n+1)*(i+0.5)*rag))))**2/(devst(i)/rag)
!print*, (ten(i)), real(n)/(log(real(primo(1)+(n+1)*(i+0.5)*rag)))
end do

  print*, stat, "stats", ((((primo(max)-primo(1))/n+1))/rag)
do i=1, ((((primo(max)-primo(1))/n+1))/rag)
  write(8,*)i, ten(i), devst(i)
end do
!p 'fort.8' u 1:2:3 w yerrorbars linestyle 1, 'fort.8' w l
!risulta perfettamente piatta

!test d'ipotesi---------------------------------------

!-------------------------------------------------------

! p 'fort.4' u 1:2 w boxes, 'fort.4' u 1:3:5 with yerrorbars linestyle 1, 'fort.4' u 1:4:5 with yerrorbars linestyle 1, 'fort.9' w l, 'fort.4' u 1:6
end program primi
