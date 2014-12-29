	SUBROUTINE EQSOLV(Ui,MEANUi,T,Ti,ums,deltat,I,RA,sum)
! input : meanui, T , Ti ,ums , sum , deltat,I,RA(K)
! output : Ui , sum

		INTEGER I,J
		DOUBLEPRECISION Ui,meanui,ums,Ti,T,deltat,sum,RA
		SUM=SUM+DEXP(I*DELTAT/Ti)*RA*deltat
		Ui=MEANUI+sqrt(2.d0*ums/Ti)*DEXP(-T/Ti)*SUM
	END SUBROUTINE EQSOLV