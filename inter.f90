program inter
  implicit none
  real, dimension(10) :: m, d, mini, maxi
  real :: s00=0, s11=0, s10=0, s01=0, s20=0, cof=0, ord=0, errm=0, errq=0, cov=0
  integer :: i, n=500
open(Unit=12, file='datimartin.txt')
  do i=1, 10
    read(12,*) m(i), d(i), mini(i), maxi(i)
    d(i)=sqrt(d(i)/n)
    write(15,*) m(i), (1./log((mini(i)+maxi(i))/2)), d(i)
  end do

do i=1, 10

s00=s00+1./(d(i)**2)
s01=s01+m(i)/(d(i)**2)
s10=s10+(1./log((mini(i)+maxi(i))/2))/(d(i)**2)
s11=s11+m(i)*(1./log((mini(i)+maxi(i))/2))/(d(i)**2)
s20=s20+(1./log((mini(i)+maxi(i))/2))**2/(d(i)**2)

end do

cof=(s00*s11-s10*s01)/(s00*s20-s10**2)
ord=(s20*s01-s10*s11)/(s00*s20-s10**2)
errm=s00/(s00*s20-s10**2)
errq=s20/(s00*s20-s10**2)
cov=-s10/(s00*s20-s10**2)

print*, "cof", cof, "ord", ord, "errm", sqrt(errm), "errq", sqrt(errq) ,"cov", cov

end program inter
