Program main
  USE DFLIB
  USE DFPORT
  implicit none
  integer(4) random1,random2,random3,random4
  integer dim_num,m,n,i,j,l,A,B,MM,mmm,nnn,ii,c,d
  integer(2) ihr, imin, isec, i100th
  integer, allocatable, dimension ( :, : ) :: element_node
  integer, allocatable, dimension ( : ) :: indx,indxy
  integer element_num
  integer element_order
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_coord,x,y,k,eps,meanu1,meanu2
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_data
  real ( kind = 8 ), allocatable, dimension ( : ) :: xx,yy,xt,kt,epst,meanu1t,meanu2t,XACOUSTIC,YACOUSTIC,rorms,dro
  real ( kind = 8 ), allocatable, dimension ( : , : , : ) :: u1,u2 !,u3 !(IF WE WANNA USE THREE DIMENTIONAL U3 IS NEEDED.)
  REAL ( KIND = 8 ), ALLOCATABLE, DIMENSION ( : , : ) ::XF,YF,EPSF,KF,MEANU1F,MEANU2F,DR,intT11,intT12,intT22,xlocal,ylocal
  integer node_data_num,p,ll,time1,last
  integer node_num
  character ( len = 80 ) :: tec_file_name = 'tecplot.txt', nameU1, nameU2
  integer tec_file_unit
  parameter stepx=1,stepy=1
  doubleprecision pi,TI,x1,x2,y1,y2,zeta1,zeta2,ums,T,deltat,R,A0,RO,sum1,sum2,W,deltax,deltay,diff11T11,diff12T12,diff22T22,retime,dv,radius,aaa,bbb,uuu,vvv
  parameter xpercent=1.0D0,ypercent=0.99D0,deltaxacoustic=0.01d0,deltayacoustic=0.01d0
  parameter tetamin=45,tetamax=45,deltateta=1
  LOGICAL(4) result1
  PARAMETER NO_OF_TIMESTEPS=1024+300 !NO_OF_TIMESTEPS must be a power of 2  + 300
  
  deltat=0.00005d0
  pi=4.d0*atan(1.d0)
  a0=340.d0
  ro=1.225d0
  
  write ( *, '(a)' )'enter the dimentions of the grid:'
  write ( *, '(a)' )'m='
  read  ( * , * )m
  write ( *, '(a)' )'n='
  read  ( * , * )n
  write (*,*) 'radius='
  read (*,*) radius
  !write ( *, '(a)' ) ' '
  !call timestamp ( )

  !write ( *, '(a)' ) ' '
  !write ( *, '(a)' ) 'TEC_IO_TEST:'
  !write ( *, '(a)' ) '  FORTRAN90 routines to read and write TEC files.'

  !write ( *, '(a)' ) ' '
  !write ( *, '(a)' ) 'TEC_IO_TEST:'
  !write ( *, '(a)' ) '  Normal end of execution.'

  !write ( *, '(a)' ) ' '
  !call timestamp ( )


  !write ( *, '(a)' ) ' '
  !write ( *, '(a)' ) 'TEC_IO_TEST01'
  !write ( *, '(a)' ) '  TEC_HEADER_READ can read the header of a TEC file.'
  !write ( *, '(a)' ) '  TEC_DATA_READ can read the data of a TEC file.'
  !write ( *, '(a)' ) ' '
  !write ( *, '(a)' ) '  In this example, we will read data from "' &
  !  // trim ( tec_file_name ) // '".'

  call tec_open_read ( tec_file_name, tec_file_unit )

  call tec_header_read ( tec_file_name, tec_file_unit, dim_num, node_num, &
    element_num, element_order, node_data_num )

  call tec_header_print ( dim_num, node_num, element_num, &
    element_order, node_data_num )

  if (m*n .ne. node_num) then
	write(*,*)'dimensions does not match.(m*n must be equal to node_num) :D'
  endif
  
  allocate ( node_coord(1:dim_num,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( node_data(1:node_data_num,1:node_num) )
  allocate ( xx(node_num),yy(m),indx(node_num),indxy(m),xt(m),kt(m),epst(m),meanu1t(m),meanu2t(m))
  allocate (x(m,n),y(m,n),k(m,n),eps(m,n),meanu1(m,n),meanu2(m,n))

  call tec_data_read ( tec_file_name, tec_file_unit, dim_num, &
    node_num, element_num, element_order, node_data_num, node_coord, &
    element_node, node_data )

  close ( unit = tec_file_unit )

  !call dmat_transpose_print_some ( dim_num, node_num, node_coord, 1, 1, &
  !  dim_num, 10, '  Coordinates of first 10 nodes:' )

  !call imat_transpose_print_some ( element_order, element_num, &
  !  element_node, 1, 1, element_order, 10, '  Nodes of first 10 elements:' )

  !call dmat_transpose_print_some ( node_data_num, node_num, node_data, 1, 1, &
  !  node_data_num, 10, '  Node data for first 10 nodes:' )

  xx(:)=node_coord(1,:)
  call SORTRX (node_num,xx,indx)
	
  do j=1,n
	do i=1,m
		yy(i)=node_coord(2,indx((j-1)*m+i))
		xt(i)=node_coord(1,indx((j-1)*m+i))
		meanu1t(i)=node_data(1,indx((j-1)*m+i))
		meanu2t(i)=node_data(2,indx((j-1)*m+i))
		kt(i)=node_data(3,indx((j-1)*m+i))
		epst(i)=node_data(4,indx((j-1)*m+i))
	end do
	call SORTRX (m,yy,indxy)
	do l=1,m
		x(l,j)=xt(indxy(l))
		y(l,j)=yy(indxy(l))
		k(l,j)=kt(indxy(l))
		eps(l,j)=epst(indxy(l))
		meanu1(l,j)=meanu1t(indxy(l))
		meanu2(l,j)=meanu2t(indxy(l))
	end do
  end do

  !--------------------
  !THIS PART REFLECTS THE QUANTITIES THAT HAVE BEEN OBTAINED FOR HALF JET TO OBTAIN THE FULL JET QUANTITIES.
  !THE QUANTITIES THAT HAVE "F" AFTER THEM ARE FOR FULL JET.
  !THE DIMENTION OF THE FULL JET IS (2*M-1,N) WE CHOOSE MM=2*M-1 SO THE DIMENTION WILL BE (MM,N)
  MM=2*M-1
  mmm=int(MM*ypercent)
  nnn=int(N*xpercent)

  ALLOCATE (XF(mmm,nnn),YF(mmm,nnn),KF(mmm,nnn),EPSF(mmm,nnn),MEANU1F(mmm,nnn),MEANU2F(mmm,nnn))
  ii=int(m*ypercent)
  DO I=1,ii
	DO J=1,nnn
		XF(I,J)=X(ii-I+1,J)
		YF(I,J)=-Y(ii-I+1,J)
		KF(I,J)=K(ii-I+1,J)
		EPSF(I,J)=EPS(ii-I+1,J)
		MEANU1F(I,J)=MEANU1(ii-I+1,J)
		MEANU2F(I,J)=-MEANU2(ii-I+1,J)
	END DO
  END DO
  DO I=1,ii+1
	DO J=1,nnn
		XF(I+ii-1,J)=X(I,J)
		YF(I+ii-1,J)=Y(I,J)
		KF(I+ii-1,J)=K(I,J)
		EPSF(I+ii-1,J)=EPS(I,J)
		MEANU1F(I+ii-1,J)=MEANU1(I,J)
		MEANU2F(I+ii-1,J)=MEANU2(I,J)		
	END DO
  END DO
  write(*,*) "Reflection of quantities was completed successfully!"
  write(*,*) 'Xmin:',xf(1,1),' Xmax:',xf(mmm,nnn)
  write(*,*) 'Ymin:',yf(1,1),' Yman:',yf(mmm,nnn)
!  call tecplot_mean_velocity_contour(Xf,Yf,meanU1f,meanU2f,mmm,nnn)  

  !--------------------
  !THIS PART MAKES THE POSITION OF THE ACOUSTIC GRID POINTS

!  MM=2*M-1
  n=tetamax-tetamin+1
  ALLOCATE (XACOUSTIC(n),YACOUSTIC(n))
  xacoustic(:)=0.0d0
  yacoustic(:)=0.0d0
  do i=1,1
	xacoustic(i)=radius*cosd(real(45))
	write(*,*) 'xacoustic:',xacoustic(i)
	yacoustic(i)=radius*sind(real(45))
	write(*,*) 'yacoustic:',yacoustic(i)
  end do
  write(*,*) "Acoustic grid is created successfully!"

  !--------------------

!-	CALL TECPLOT_VELOCITY_COTOUR(XF,YF,meanU1F,meanU2F,MM,N,NO_OF_TIMESTEPS)	

  deallocate ( node_coord,indx,indxy )
  deallocate ( element_node,xx,yy )
  deallocate ( node_data,xt,kt,epst,meanu1t,meanu2t )
  DEALLOCATE ( X,Y,K,EPS,MEANU1,MEANU2)
  
  !--------------------
  ! THIS PART GENERATES THE VELOCITY FLUCTUATIONS FROM THE MEAN FLOW QUANTITIES USING THE G.AHMADI METHOD
  ! FOR GENERATING TOW UNIFORM RANDOM NUMBER, WE USE THE COMPUTR CLOCK SO EVERY TIME DIFFERENT RANDOM NUMBERS WILL BE GENERATED
  ! GUSSIAN RANDOM NUMBERS WILL BE GENERATED IN THE SUBRUTINE GAUSSIAN() BY USING THESE TWO UNIFORM RANDOM NUMBERS
	
	result1 = MAKEDIRQQ('e:\DATA')
	result1 = CHANGEDIRQQ ('e:\data')
	do j=1,nnn
		write (nameu1 ,*) j
		result1 = MAKEDIRQQ(nameu1)
!		IF (result1) THEN
!		   WRITE (*,*) 'New subdirectory successfully created',j
!		ELSE
!		   WRITE (*,*) 'Failed to create subdirectory'
!		END IF
	end do

	meanu2f(:,:)=0.d0
	random3=98765
	random4=56789
	do j=1,nnn
		do i=1,mmm
			call GETTIM (ihr, imin, isec, i100th)
			if (mod(i100th,2).eq.0) then
				aaa=ran(random3)
				bbb=ran(random4)
				random1=(random3/100)*(i100th+1)
				random2=(random4/100)*(i100th+1)
			else
				aaa=ran(random3)
				bbb=ran(random4)
				random1=(random3/100)*i100th
				random2=(random4/100)*i100th
			end if
			t=0.0d0

			SUM1=0.D0
			SUM2=0.D0
			write (nameU1 ,*) 'E:\DATA\',j,'\U1_',i,'_',j,'.dat'
!			write (*,*) nameu1
			open (unit= 111, file=nameU1, form = 'binary')
			write (nameU2 ,*) 'E:\DATA\',j,'\U2_',i,'_',j,'.dat'
!			write (*,*) nameu2
			open (unit= 112, file=nameU2, form = 'binary')
			do l=1,NO_OF_TIMESTEPS
				T=T+deltat

1				X1=2.D0*RAN(RANDOM1)-1.D0
				X2=2.D0*RAN(RANDOM2)-1.D0
				W=X1*X1+X2*X2
				IF (W.GE.1.0) GOTO 1
				W=SQRT((-2.0D0*LOG(W))/W)
				Y1=X1*W
				Y2=X2*W

				call TTI(kF(int(mmm/2),1),epsF(int(mmm/2),1),TI)
				call ZETAi(y1,deltat,zeta1)
				call ZETAi(y2,deltat,zeta2)
				call UiMS(kF(i,j),ums) !kF(int(mmm/2),1)
				call EQSOLV(UUU,meanu1F(I,J),T,TI,ums,deltat,L,ZETA1,SUM1)
				call EQSOLV(VVV,meanu2F(I,J),T,TI,ums,deltat,L,ZETA2,SUM2)
				write (111) UUU
				write (112) VVV
			end do
			close (111)
			close (112)
		end do
	end do

	write(*,*) "Velocity fluctuations was generated successfully!"
!-----------------------------------------------

!	CALL TECPLOT_VELOCITY_COTOUR(XF,YF,U1,U2,mmm,nnn,NO_OF_TIMESTEPS)
!	CALL TECPLOT_VELOCITY_FLUCTUATION_COTOUR(XF,YF,U1,U2,mmm,nnn,NO_OF_TIMESTEPS,MEANU1F,MEANU2F)
!	write(*,*) 'velocity data was written successfully!'

	DEALLOCATE ( KF,EPSF,MEANU1F,MEANU2F )
	allocate (DR(n,no_of_timesteps),intT11(5,5),intT12(5,5),intT22(5,5),dro(no_of_timesteps))
	allocate (xlocal(5,5),ylocal(5,5))
    allocate (u1(mmm,nnn,300),u2(mmm,nnn,300)) !,u3(mmm,nnn,NO_OF_TIMESTEPS))
!-----------------------------------------------
	
  open(unit=2,file='timehistoryofro.dat')
   dr(:,:)=0.0d0
   dro(:)=0.0d0
	last=300  
	do j=1,nnn
		do i=1,mmm
			write (nameU1 ,*) 'E:\DATA\',j,'\U1_',i,'_',j,'.dat'
			open (unit= 111, file=nameU1, form = 'binary')
			write (nameU2 ,*) 'E:\DATA\',j,'\U2_',i,'_',j,'.dat'
			open (unit= 112, file=nameU2, form = 'binary')
			
			do time1 = 1,last
				read (111) u1(i,j,time1)
				read (112) u2(i,j,time1)
			end do
			close (111)
			close (112)
		end do
	end do
	
	P=0
   do l=1,no_of_timesteps-last

 	If (l.gt.Last*(p+1)) then
	do j=1,nnn
		do i=1,mmm
			write (nameU1 ,*) 'E:\DATA\',j,'\U1_',i,'_',j,'.dat'
			open (unit= 111, file=nameU1, form = 'binary')
			write (nameU2 ,*) 'E:\DATA\',j,'\U2_',i,'_',j,'.dat'
			open (unit= 112, file=nameU2, form = 'binary')

			
			do time1 = 1,(P+1)*250
				read (111) uuu
				read (112) vvv
			end do
			do time1 = 1,last
				read (111) u1(i,j,time1)
				read (112) u2(i,j,time1)
			end do
			close (111)
			close (112)
		end do
	end do
	
	P=P+1
	end if

	ll=l-P*300
	do i=1,1,1 !tetamin,tetamax,deltateta
	intT11(:,:)=0.d0
	intT12(:,:)=0.d0
	intT22(:,:)=0.d0	  
     do c=1,5
	  do d=1,5
	   xlocal(c,d)=xacoustic(i)+(real(d-3)*deltaxacoustic)
	   ylocal(c,d)=yacoustic(i)+(real(c-3)*deltayacoustic)
	  end do
	 end do

	 do c=1,5
	  do d=1,5
	   x1=xlocal(c,d)
	   x2=ylocal(c,d)
	   do a=1,nnn,stepx
	    do b=1,mmm,stepy
		 Y1=XF(b,a)
		 Y2=YF(b,a)
		 if ((b+stepy).lt.mmm) then
			deltay=abs(yf(b+stepy,a)-yf(b,a))
		 else
			deltay=abs(yf(b,a)-yf(b-stepy,a))
		 end if
		 if ((a+stepx).lt.nnn) then
			deltax=abs(xf(b,a+stepx)-xf(b,a))
		 else
			deltax=abs(xf(b,a)-xf(b,a-stepx))
		 end if
		 dv=deltax*deltay
	     r=sqrt((y1-x1)*(y1-x1)+(y2-x2)*(y2-x2))
		 retime=ll-INT(r/a0/deltat)               
		 if ((r.gt.0).and.(retime.gt.0)) then
		  intT11(c,d)=intT11(c,d)+(ro*u1(b,a,retime)*u1(b,a,retime))/r*dv
		  intT12(c,d)=intT12(c,d)+(ro*u1(b,a,retime)*u2(b,a,retime))/r*dv
		  intT22(c,d)=intT22(c,d)+(ro*u2(b,a,retime)*u2(b,a,retime))/r*dv
		 end if
		end do
	   end do
   	  end do
	 end do

	diff11T11=(-intT11(3,1)+16.d0*intT11(3,2)-30.d0*intT11(3,3)+16.d0*intT11(3,4)-intT11(3,5))/12.d0/deltaxacoustic**2.0d0
	diff22T22=(-intT11(1,3)+16.d0*intT11(2,3)-30.d0*intT11(3,3)+16.d0*intT11(4,3)-intT11(5,3))/12.d0/deltayacoustic**2.0d0
	diff12T12=((intT12(1,1)-8.d0*intT12(2,1)+8.d0*intT12(4,1)-intT12(5,1))	&
		 -8.d0*(intT12(1,2)-8.d0*intT12(2,2)+8.d0*intT12(4,2)-intT12(5,2))	&
		 +8.d0*(intT12(1,4)-8.d0*intT12(2,4)+8.d0*intT12(4,4)-intT12(5,4))	&
		 -1.0d0*(intT12(1,5)-8.d0*intT12(2,5)+8.d0*intT12(4,5)-intT12(5,5)))	&
		 /144.d0/deltaxacoustic/deltayacoustic
	dr(i,l)=(diff11T11+2.d0*diff12T12+diff22T22)/4.0d0/pi/a0**2.d0
	write(*,*) dr(i,l)
	write(2,7) dr(i,l)
	dro(l)=dr(i,l)
   end do
!   write(*,*) real(l)*real(deltat),'/',real(no_of_timesteps)*real(deltat),l,'/',no_of_timesteps
  end do
7 format(e20.10)
close(2)

!----------------------------------------------
!call FFT SUBROUTINE
!----------------------------------------------

call realft(dro,no_of_timesteps-last,1)
open(unit=3,file='FFTtimehistoryofro.dat')
do l=50,562
 write(3,3) l,abs(dro(l)) 
end do
close(3)
3 format(I8,1x,e20.10)



!----------------------------------------------
!allocate(rorms(n))
!rorms(:)=0.d0
!do i=tetamin,tetamax,deltateta
! do l=abs(no_of_timesteps/2),no_of_timesteps
!  rorms(i)=rorms(i)+dr(i,l)**2
! end do
! rorms(i)=sqrt(rorms(i)/no_of_timesteps)
!end do
!open(unit=1,file='rorms.dat')
!do i=tetamin,tetamax,deltateta
! write(1,*) i,rorms(i)
!end do
!close(1)

!-----------------------------------------------
!call beeps
stop
end

!-----------------------------------------------

!subroutine tecplot(dr,no_of_timesteps,n,deltateta)

!integer i,n,no_of_timesteps,deltateta
!doubleprecision dr(n,no_of_timesteps)

!open(unit=7,file='density_fluctuation.dat')
!write(7,*) 'variables=x,','"','dr','"'
!do l=1,no_of_timesteps
!	write(7,*) ' Zone T =', '"' ,'t',l, '"', 'I =', int(n/2)
!	do i=1,n,deltateta
!		write(7,2) i,dr(i,l)
!	end do
!end do
!close(7)
!2 format (i5,2x,e20.10)

!end subroutine tecplot

!-----------------------------------------------

!SUBROUTINE TECPLOT_VELOCITY_COTOUR(X,Y,U1,U2,M,N,L)
! THIS SUBROUTINE SAVES A FILE FOR TECPLOT
!INTEGER M,N,L,I,J
!DOUBLEPRECISION X(M,N),Y(M,N),U1(M,N,L),U2(M,N,L) !,U1(M,N),U2(M,N) 

!	OPEN(UNIT=5,FILE='VELOCITY.DAT')
!	WRITE(5,*) 'TITLE = "VELOCIY CONTOURS"'
!	WRITE(5,*) 'variables = "x", "y", "U1", "U2", "U"'
!	WRITE(5,*) 'zone i=',M,' , j= ',N,' ,DATAPACKING=POINT'
!	WRITE(*,*) 'i=',M,' ,j= ',N
!	DO J=1,N
!		DO I=1,M
!			 WRITE(5,5) REAL(X(I,J)),REAL(Y(I,J)),REAL(U1(I,J,L)),REAL(U2(I,J,L)),REAL(SQRT(U1(I,J,L)**2+U2(I,J,L)**2)) !,REAL(U1(I,J)),REAL(U2(I,J)),REAL(SQRT(U1(I,J)**2+U2(I,J)**2)) !
!		END DO
!	END DO
!5	FORMAT(2X,d20.10,2x,d20.10,2x,d20.10,2x,d20.10,2x,d20.10,2x)
!	CLOSE(5)
	
!END SUBROUTINE TECPLOT_VELOCITY_COTOUR

!-----------------------------------------------

!SUBROUTINE TECPLOT_VELOCITY_FLUCTUATION_COTOUR(X,Y,U1,U2,M,N,L,MEANU1F,MEANU2F)
! THIS SUBROUTINE SAVES A FILE FOR TECPLOT
!INTEGER M,N,L,I,J
!DOUBLEPRECISION X(M,N),Y(M,N),U1(M,N,L),U2(M,N,L),MEANU1F(M,N),MEANU2F(M,N) !,U1(M,N),U2(M,N) 

!	OPEN(UNIT=5,FILE='VELOCITY_fluctuations.DAT')
!	WRITE(5,*) 'TITLE = "VELOCIY CONTOURS"'
!	WRITE(5,*) 'variables = "x", "y", "U1", "U2", "U"'
!	WRITE(5,*) 'zone i=',M,' , j= ',N,' ,DATAPACKING=POINT'
!	WRITE(*,*) 'i=',M,' ,j= ',N
!	DO J=1,N
!		DO I=1,M
!			 WRITE(5,5) REAL(X(I,J)),REAL(Y(I,J)),REAL((U1(I,J,L)-MEANU1F(I,J))),REAL((U2(I,J,L)-MEANU2F(I,J))),REAL((SQRT((U1(I,J,L)-MEANU1F(I,J))**2+(U2(I,J,L)-MEANU2F(I,J))**2)))
!		END DO
!	END DO
!5	FORMAT(2X,d20.10,2x,d20.10,2x,d20.10,2x,d20.10,2x,d20.10,2x)
!	CLOSE(5)
	
!END SUBROUTINE TECPLOT_VELOCITY_FLUCTUATION_COTOUR

!-----------------------------------------------

!SUBROUTINE TECPLOT_mean_VELOCITY_CONTOUR(X,Y,U1,U2,M,N)
! THIS SUBROUTINE SAVES A FILE FOR TECPLOT
!INTEGER M,N,I,J
!DOUBLEPRECISION X(M,N),Y(M,N),U1(M,N),U2(M,N) 

!	OPEN(UNIT=5,FILE='Mean_VELOCITY.DAT')
!	WRITE(5,*) 'TITLE = "MEAN VELOCIY CONTOURS"'
!	WRITE(5,*) 'variables = "x", "y", "U1", "U2", "U"'
!	WRITE(5,*) 'zone i=',M,' , j= ',N,' ,DATAPACKING=POINT'
!	WRITE(*,*) 'i=',M,' ,j= ',N
!	DO J=1,N
!		DO I=1,M
!			 WRITE(5,5) REAL(X(I,J)),REAL(Y(I,J)),REAL(U1(I,J)),REAL(U2(I,J)),REAL(SQRT(U1(I,J)**2+U2(I,J)**2))
!		END DO
!	END DO
!5	FORMAT(2X,d20.10,2x,d20.10,2x,d20.10,2x,d20.10,2x,d20.10,2x)
!	CLOSE(5)
	
!END SUBROUTINE TECPLOT_mean_VELOCITY_CONTOUR

!-----------------------------------------------

!Subroutine beeps

!integer i
!do i=600,300,-20
!call beepqq(i,100)
!end do
!do i=320,600,20
!call beepqq(i,100)
!end do

!end subroutine beeps