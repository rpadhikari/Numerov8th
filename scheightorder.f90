program schrodinger
use params
  call allocarrays

  open(1,file='potential.dat',action='write')
  do i=1, mesh
    x(i)=xmin+i*dx
    v(i)=De*(1.0d0-(exp((xe-x(i))/a)-1.0d0)**2)
    !v(i)=De*(exp(-2.0d0*a*(x(i)-xe))-2.0d0*(exp(-a*(x(i)-xe))))
    write(1,*) x(i), v(i)
  end do
  close(1)


  do i=1, mesh
    do j=1, mesh
    !---------------------------------------------------------------------
    !Declaring Vmatrix, Amatrix and Bmatrix
      if(i==j) then
        Vmat(i,i)=v(i)
        Amat(i,i)=-6.0d0
        Bmat(i,i)=166768.0d0
      elseif(i==j+1 .or. i==j-1) then
        Amat(i,j)=1.0d0
        Bmat(i,j)=178108.0d0
      elseif(i==j+2 .or. i==j-2) then
        Amat(i,j)=1.0d0
        Bmat(i,j)=7112.0d0
      elseif(i==j-3 .or. i==j+3) then
        Amat(i,j)=1.0d0
        Bmat(i,j)=13636.0d0
      else
        Vmat(i,j)=0.0d0
        Amat(i,j)=0.0d0
        Bmat(i,j)=0.0d0
      end if

    end do
  enddo
!---------------------------------------------------------------------------
  Bmat=Bmat/40320.0d0
  Amat=-Cse*Amat
  Amat=Amat/dx2

  BmatInv=Bmat

  call DGETRF(M, N, BmatInv, LDA, IPIV, INFO)

  call DGETRI(N, BmatInv, LDA, IPIV, WORK1, LWORK1, INFO)

  Tmat=matmul(BmatInv, Amat)

  Hmat=Tmat+Vmat

!----------------------------------------------------------------------------
  call DSYEV(JOBZ, UPLO, N, Hmat, LDA, W, WORK2, LWORK2, INFO)

  Psi=Hmat

  evalcount=0
  open(2,file='eigval.dat',action='write')
  open(3,file='psi.dat',action='write')
  open(4,file='prob.dat',action='write')
  do i=1, mesh
    if (W(i) .gt. De .and. W(i) .lt. 0) then
      write(2,*) i, W(i)
      isecond=i
      evalcount=evalcount+1
      normF=dot_product(Psi(:,i),Psi(:,i))*dx
      Psi(:,i)=Psi(:,i)/sqrt(normF)
      Prob(:,i)=Psi(:,i)*Psi(:,i)
    end if
  end do
  close(2)

  open(2,file='eigval.dat',action='read')
  read(2,*) ifirst
  close(2)

  do j=1, mesh
    write(3,*) x(j), Psi(j,ifirst:isecond)
    write(4,*) x(j), Prob(j,ifirst:isecond)
  end do
  close(3)
  close(4)

  call deallocarrays
end program schrodinger
