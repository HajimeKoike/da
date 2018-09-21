!==============================================================!
! multi model Serial Ensemble Square Root Filter on L63 model
! [Reference]
!   Lorenz, 1963
!   Whitaker and Hamill, 2002
!   Otsuka and Miyoshi, 2015
! [History]
!   Hajime KOIKE created,  Aug. 31st, 2018
!   Hajime KOIKE modified, Sep. 6th, 2018
!==============================================================!

PROGRAM Lorenz63_multiSeSRF_DBF
  IMPLICIT NONE

  !-----------------------------------------
  ! PARAMETERS
  !-----------------------------------------
  INTEGER,PARAMETER :: nt_asm=1000, nt_prd=1000
  INTEGER,PARAMETER :: obs_interval=100
  REAL(8),PARAMETER :: dt=1.0d-3

  REAL(8),PARAMETER :: PI=3.14159265358979d0

  INTEGER,PARAMETER :: q=3 !model dimension DO NOT CHANGE
  INTEGER,PARAMETER :: n=5 !models

  INTEGER,PARAMETER :: Nall=10000 !total ensemble size
  !SUM(m(0:n,t/oobs_interval)=Nall
  INTEGER,PARAMETER :: p=3 !observation dimension

  !-----------------------------------------
  ! STATE VARIABLES
  !-----------------------------------------
  REAL(8) :: x_t(0:nt_asm+nt_prd)
  REAL(8) :: y_t(0:nt_asm+nt_prd)
  REAL(8) :: z_t(0:nt_asm+nt_prd)

  REAL(8) :: x_s(0:nt_asm+nt_prd,n)
  REAL(8) :: y_s(0:nt_asm+nt_prd,n)
  REAL(8) :: z_s(0:nt_asm+nt_prd,n)

  REAL(8) :: parm2(n)
  INTEGER :: m(0:nt_asm,n)
  !------------ensemble members------------
  !
  ! x_da_k(t,:) = x_da(t) + x_prtb(t,:)
  !
  !----------------------------------------
  REAL(8),ALLOCATABLE :: x_da_k(:,:)
  REAL(8),ALLOCATABLE :: y_da_k(:,:)
  REAL(8),ALLOCATABLE :: z_da_k(:,:)
!
!  !------------ensemble mean---------------
  REAL(8) :: x_da(0:nt_asm+nt_prd,n)
  REAL(8) :: y_da(0:nt_asm+nt_prd,n)
  REAL(8) :: z_da(0:nt_asm+nt_prd,n)

  !----------------------------------------
  ! WEIGHTS, LIKELIHOODS
  !----------------------------------------
  REAL(8) :: W(0:nt_asm/obs_interval,n)
  REAL(8) :: Pr(0:nt_asm/obs_interval,n)

  !-----------------------------------------
  ! OBSERVATION VARIABLES
  !-----------------------------------------
  REAL(8) :: Y(0:nt_asm/obs_interval,p)

  !-----------------------------------------
  ! OUTPUT CONTROL
  !-----------------------------------------
  INTEGER,PARAMETER :: output_interval=50
  CHARACTER(LEN=7) :: obs_chr(0:nt_asm)

  !-----------------------------------------
  ! MATRIX
  !-----------------------------------------
  REAL(8) :: Ef(q,Nall)
  REAL(8) :: Pf(q,q)
  REAL(8) :: Pa(q,q)
  REAL(8) :: R(p,p)
  REAL(8) :: Kg(q,1)
  REAL(8) :: H(p,q)
  REAL(8) :: Id(q,q) ! identity matrix

  !-----------------------------------------
  ! WORKING VARIABLES
  !-----------------------------------------
  INTEGER :: i,j,k,k4,l,s,t
  REAL(8) :: X_mean(q,1),dX(q,Nall)
  REAL(8) :: denom(1,1)
  REAL(8) :: alpha(1,1)
  REAL(8) :: innov(1,1)
  REAL(8) :: noise1,noise2
  REAL(8) :: gnoise
  REAL(8) :: H8(1,q)
  REAL(8) :: PfHt(q,1)
  REAL(8) :: Hx(p,n)
  REAL(8) :: cpr

  !-----------------------------------------
  ! RMSEs
  !-----------------------------------------
  REAL(8) :: SQ_s(-1:nt_asm+nt_prd,q),SQ_da(-1:nt_asm+nt_prd,q)
  REAL(8) :: RMSE_s(0:nt_asm+nt_prd,q),RMSE_da(0:nt_asm+nt_prd,q)


  CALL RANDOM_SEED()

  !-----------------------------------------
  ! INIT VALUE SETTING
  ! x_t,y_t,z_t
  ! x_s,y_s,z_s
  ! Pa
  !-----------------------------------------

  x_t(0)=10.0d0;y_t(0)=14.0d0;z_t(0)=24.0d0
  x_s(0,1:n)=11.0d0;y_s(0,1:n)=13.0d0;z_s(0,1:n)=25.0d0

  Pf(1,1)=1.0d0;Pf(1,2)=0.0d0;Pf(1,3)=0.0d0
  Pf(2,1)=0.0d0;Pf(2,2)=1.0d0;Pf(2,3)=0.0d0
  Pf(3,1)=0.0d0;Pf(3,2)=0.0d0;Pf(3,3)=1.0d0

  Pa=Pf

  parm2(1)=27.0d0
  parm2(2)=27.5
  parm2(3)=28.0d0
  parm2(4)=28.5
  parm2(5)=29.0d0

  !-----------------------------------------
  ! DEFINE R,H
  !-----------------------------------------
  R(1,1)=0.1d0;R(1,2)=0.0d0;R(1,3)=0.0d0
  R(2,1)=0.0d0;R(2,2)=0.1d0;R(2,3)=0.0d0
  R(3,1)=0.0d0;R(3,2)=0.0d0;R(3,3)=0.1d0

  H(1,1)=1.0d0;H(1,2)=0.0d0;H(1,3)=0.0d0
  H(2,1)=0.0d0;H(2,2)=1.0d0;H(2,3)=0.0d0
  H(3,1)=0.0d0;H(3,2)=0.0d0;H(3,3)=1.0d0


  DO i=1,q
    DO j=1,q
      IF (i .EQ. j) THEN
        Id(i,j)=1.0d0
      ELSE
        Id(i,j)=0.0d0
      ENDIF
    ENDDO
  ENDDO
  !-----------------------------------------
  ! NATURE RUN, OBSERVATION
  !-----------------------------------------

  DO t=1,nt_asm+nt_prd
    CALL tinteg_L63_true(x_t(t-1),y_t(t-1),z_t(t-1),&
      &x_t(t),y_t(t),z_t(t))

    IF (MOD(t, obs_interval) == 0 .AND. (t <= nt_asm)) THEN
      CALL RANDOM_NUMBER(noise1)
      CALL RANDOM_NUMBER(noise2)
      gnoise=SQRT(R(1,1))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
        &*COS(2.0d0*PI*noise2)
      Y(t/obs_interval,1)=x_t(t)+gnoise

      CALL RANDOM_NUMBER(noise1)
      CALL RANDOM_NUMBER(noise2)
      gnoise=SQRT(R(2,2))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
        &*COS(2.0d0*PI*noise2)
      Y(t/obs_interval,2)=y_t(t)+gnoise

      CALL RANDOM_NUMBER(noise1)
      CALL RANDOM_NUMBER(noise2)
      gnoise=SQRT(R(3,3))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
        &*COS(2.0d0*PI*noise2)
      Y(t/obs_interval,3)=z_t(t)+gnoise

    ENDIF
  ENDDO

  !-----------------------------------------
  ! SIMULATION WITHOUT DA
  !-----------------------------------------

  DO t=1,nt_asm+nt_prd
    DO s=1,n
      CALL tinteg_L63_sim(x_s(t-1,s),y_s(t-1,s),z_s(t-1,s),&
        &x_s(t,s),y_s(t,s),z_s(t,s),parm1=10.0d0,parm2=parm2(s),&
        &parm3=8.0d0/3.0d0)
    ENDDO
  ENDDO

  !-----------------------------------------
  ! SIMULATION WITH DA
  !-----------------------------------------

  DO s=1,n
    x_da(0,s)=x_s(0,s)
    y_da(0,s)=y_s(0,s)
    z_da(0,s)=z_s(0,s)
  ENDDO

  m(0,:)=Nall/n

  ALLOCATE(x_da_k(Nall,n),y_da_k(Nall,n),z_da_k(Nall,n))

  DO s=1,n
    DO k=1,m(0,s)
      CALL RANDOM_NUMBER(noise1)
      CALL RANDOM_NUMBER(noise2)
      gnoise=SQRT(Pf(1,1))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
        &*COS(2.0d0*PI*noise2)
      x_da_k(k,s)=x_da(0,s)+gnoise

      CALL RANDOM_NUMBER(noise1)
      CALL RANDOM_NUMBER(noise2)
      gnoise=SQRT(Pf(2,2))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
        &*COS(2.0d0*PI*noise2)
      y_da_k(k,s)=y_da(0,s)+gnoise

      CALL RANDOM_NUMBER(noise1)
      CALL RANDOM_NUMBER(noise2)
      gnoise=SQRT(Pf(3,3))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
        &*COS(2.0d0*PI*noise2)
      z_da_k(k,s)=z_da(0,s)+gnoise
    ENDDO
  ENDDO

  !-----------------------------------------------------------
  ! DA: EnSRF state update & DBF emsemble size optimization
  !-----------------------------------------------------------
  Pr(0/obs_interval,:)=1.0d0/n

  DO t=1,nt_asm
    DO s=1,n
      DO k=1,m(t-1,s)
          CALL tinteg_L63_sim(x_da_k(k,s),y_da_k(k,s),z_da_k(k,s),x_da_k(k,s),y_da_k(k,s),z_da_k(k,s),&
            &parm1=10.0d0,parm2=parm2(s),parm3=8.0d0/3.0d0)
      ENDDO
    ENDDO

    WRITE(*,*) ">>>FORECAST"
    DO s=1,n
      x_da(t,s)=SUM(x_da_k(1:m(t-1,s),s))/m(t-1,1)
      y_da(t,s)=SUM(y_da_k(1:m(t-1,s),s))/m(t-1,1)
      z_da(t,s)=SUM(z_da_k(1:m(t-1,s),s))/m(t-1,1)
      WRITE(*,'(3F12.7)') x_da(t,s),y_da(t,s),z_da(t,s)
    ENDDO

    IF (MOD(t, obs_interval) == 0) THEN

      X_mean(1,1)=SUM(x_da(t,1:n))/n
      X_mean(2,1)=SUM(y_da(t,1:n))/n
      X_mean(3,1)=SUM(z_da(t,1:n))/n

      k4=0
      DO s=1,n
        DO k=1,m(t-1,s)
          k4=k4+1
          dX(1,k4)=X_mean(1,1)-x_da_k(k,s)
          dX(2,k4)=X_mean(2,1)-y_da_k(k,s)
          dX(3,k4)=X_mean(3,1)-z_da_k(k,s)
        ENDDO
      ENDDO

!      WRITE(*,*) "----shape of dX should be (n,Nall) = (3,100)----"
!      WRITE(*,'(2I5)') SHAPE(dX)

      Ef = SQRT(1.0d0/(Nall-1))*dX
      Pf=MATMUL(Ef,TRANSPOSE(Ef))
      WRITE(*,*)
      WRITE(*,*) "-----Pf-----"
      WRITE(*,'(3F10.7)') Pf

      !--------------------------------------
      ! DA: SeSRF
      !--------------------------------------
      DO l=1,p
        H8(1,1:q)=H(l,1:q)
        PfHt=MATMUL(Pf,TRANSPOSE(H8))
        denom=1.0d0/(MATMUL(H8,PfHt)+R(l,l))
        WRITE(*,*) "----denom----"
        WRITE(*,'(F7.2)') denom
        Kg=MATMUL(PfHt,denom)
        WRITE(*,*) "----Kalman gain----"
        WRITE(*,'(F7.2)') Kg
        alpha=1.0d0/(1.0d0+SQRT(R(l,l)*denom))
        WRITE(*,*) "----alpha----"
        WRITE(*,'(F7.2)') alpha
        innov=Y(t/obs_interval,l)-MATMUL(H8,X_mean)

        X_mean=X_mean+MATMUL(Kg,innov)

        DO k=1,Nall
          dX(1:q,k)=MATMUL((Id-alpha(1,1)*MATMUL(Kg,H8)),&
            &dX(1:q,k))
        ENDDO
      ENDDO

      k4=0
      DO s=1,n
        DO k=1,m(t-1,s)
          k4=k4+1
          x_da_k(k,s)=X_mean(1,1)-dX(1,k4)
          y_da_k(k,s)=X_mean(2,1)-dX(2,k4)
          z_da_k(k,s)=X_mean(3,1)-dX(3,k4)
        ENDDO
      ENDDO
      WRITE(*,*) ">>>ANALYSIS"
      DO s=1,n
        x_da(t,s)=SUM(x_da_k(1:m(t-1,s),s))/m(t-1,s)
        y_da(t,s)=SUM(y_da_k(1:m(t-1,s),s))/m(t-1,s)
        z_da(t,s)=SUM(z_da_k(1:m(t-1,s),s))/m(t-1,s)
        WRITE(*,'(3F12.7)') x_da(t,s),y_da(t,s),z_da(t,s)
      ENDDO
      !----------------------------------------
      ! DA:Discrete Bayesian Filter
      !----------------------------------------
      DO s=1,n
        Hx(:,s)=0.0d0
        DO i=1,p
          Hx(i,s)=Hx(i,s)+H(i,1)*x_da(t,s)+H(i,2)*y_da(t,s)+H(i,3)*z_da(t,s)
        ENDDO
        W(t/obs_interval,s)=1.0d0/SQRT(SUM((Hx(:,s)-Y(t/obs_interval,:))**2)/p)
      ENDDO
      WRITE(*,*)
      WRITE(*,*) "---------Likelihood---------"
      WRITE(*,'(5F25.12)') W(t/obs_interval,:)
      DO s=1,n
        Pr(t/obs_interval,s)=W(t/obs_interval,s)*Pr(t/obs_interval-1,s)
      ENDDO
        Pr(t/obs_interval,:)=Pr(t/obs_interval,:)/SUM(Pr(t/obs_interval,1:n))
      WRITE(*,*)
      WRITE(*,*) "---------Probability---------"
      WRITE(*,'(5F25.12)') Pr(t/obs_interval,:)

      m(t,:)=0
      DO k4=1,Nall
        CALL RANDOM_NUMBER(noise1)
        cpr=0.0d0
        s=0
        DO WHILE (cpr .LE. noise1)
          s=s+1
          cpr=cpr+Pr(t/obs_interval,s)
        ENDDO
        m(t,s)=m(t,s)+1
      ENDDO
      WRITE(*,*)
      WRITE(*,*) "-------optimal ensemble sizes--------"
      WRITE(*,'(5I5)') m(t,:)

      DEALLOCATE(x_da_k,y_da_k,z_da_k)
      ALLOCATE(x_da_k(Nall,n),y_da_k(Nall,n),z_da_k(Nall,n))
      DO s=1,n
        DO k=1,m(t,s)
          CALL RANDOM_NUMBER(noise1)
          CALL RANDOM_NUMBER(noise2)
          gnoise=SQRT(Pf(1,1))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
            &*COS(2.0d0*PI*noise2)
          x_da_k(k,s)=x_da(t,s)+gnoise

          CALL RANDOM_NUMBER(noise1)
          CALL RANDOM_NUMBER(noise2)
          gnoise=SQRT(Pf(2,2))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
            &*COS(2.0d0*PI*noise2)
          y_da_k(k,s)=y_da(t,s)+gnoise

          CALL RANDOM_NUMBER(noise1)
          CALL RANDOM_NUMBER(noise2)
          gnoise=SQRT(Pf(3,3))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
            &*COS(2.0d0*PI*noise2)
          z_da_k(k,s)=z_da(t,s)+gnoise
        ENDDO
      ENDDO
    ELSE
      DO s=1,n
        m(t,s)=m(t-1,s)
      ENDDO
    ENDIF
  ENDDO

  DO t=nt_asm+1,nt_asm+nt_prd
    DO s=1,n
      CALL tinteg_L63_sim(x_da(t-1,s),y_da(t-1,s),z_da(t-1,s),x_da(t,s),y_da(t,s),z_da(t,s),&
        &parm1=10.0d0,parm2=parm2(s),parm3=8.0d0/3.0d0)
    ENDDO
  ENDDO

  !-----------------------------------------
  ! OUTPUT
  !-----------------------------------------

  obs_chr(0:nt_asm)="No obs"
  DO t=1,nt_asm
    IF (MOD(t,obs_interval)==0) THEN
      WRITE(obs_chr(t),'(F7.3)') Y(t/obs_interval,1)
    ENDIF
  ENDDO
  WRITE(*,*)
  WRITE(*,*) "================================================"
  WRITE(*,*) " multimodel SeSRF with Discrete Bayesian Filter "
  WRITE(*,*) "================================================"
  WRITE(*,*)
  WRITE(*,'(A,F7.2,A,F7.2)') "Assimilation Period: t=", &
    & 0.0,"-",dt*nt_asm
  WRITE(*,'(A,F7.2,A,F7.2)') "Prediction Period: t=", &
    & dt*nt_asm,"-",dt*(nt_asm+nt_prd)
  WRITE(*,*)
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "        Assimilation Period:x         "
  WRITE(*,*) "  [Time]   [True]   [No DA 2] [No DA 3] [No DA 4]  [DA 2]   [DA 3]   [DA 4]   [Obs]"
  DO t=0,nt_asm
    IF (MOD(t,output_interval)==0) THEN
      WRITE(*,'(F7.2,7F10.3,4X,A)') dt*t,x_t(t),x_s(t,2),x_s(t,3),x_s(t,4),&
        &x_da(t,2),x_da(t,3),x_da(t,4),obs_chr(t)
    ENDIF
  ENDDO
  DO t=nt_asm+1,nt_asm+nt_prd
    IF(MOD(t,output_interval)==0) THEN
      WRITE(*,'(F7.2,7F10.3)') dt*t,x_t(t),x_s(t,2),x_s(t,3),x_s(t,4),&
        &x_da(t,2),x_da(t,3),x_da(t,4)
    ENDIF
  ENDDO

  obs_chr(0:nt_asm)="No obs"
  DO t=1,nt_asm
    IF (MOD(t,obs_interval)==0) THEN
      WRITE(obs_chr(t),'(F7.3)') Y(t/obs_interval,2)
    ENDIF
  ENDDO
  WRITE(*,*)
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "        Assimilation Period:y        "
  WRITE(*,*) "  [Time]   [True]   [No DA 2] [No DA 3] [No DA 4]  [DA 2]   [DA 3]   [DA 4]   [Obs]"
  DO t=0,nt_asm
    IF (MOD(t,output_interval)==0) THEN
      WRITE(*,'(F7.2,7F10.3,4X,A)') dt*t,y_t(t),y_s(t,2),y_s(t,3),y_s(t,4),&
        &y_da(t,2),y_da(t,3),y_da(t,4),obs_chr(t)
    ENDIF
  ENDDO
  DO t=nt_asm+1,nt_asm+nt_prd
    IF(MOD(t,output_interval)==0) THEN
      WRITE(*,'(F7.2,7F10.3)') dt*t,y_t(t),y_s(t,2),y_s(t,3),y_s(t,4),&
        &y_da(t,2),y_da(t,3),y_da(t,4)
    ENDIF
  ENDDO

  obs_chr(0:nt_asm)="No obs"
  DO t=1,nt_asm
    IF (MOD(t,obs_interval)==0) THEN
      WRITE(obs_chr(t),'(F7.3)') Y(t/obs_interval,3)
    ENDIF
  ENDDO
  WRITE(*,*)
  WRITE(*,*) "--------------------------------------"
  WRITE(*,*) "        Assimilation Period:z         "
  WRITE(*,*) "  [Time]   [True]   [No DA 2] [No DA 3] [No DA 4]  [DA 2]   [DA 3]   [DA 4]   [Obs]"
  DO t=0,nt_asm
    IF (MOD(t,output_interval)==0) THEN
      WRITE(*,'(F7.2,7F10.3,4X,A)') dt*t,z_t(t),z_s(t,2),z_s(t,3),z_s(t,4),&
        &z_da(t,2),z_da(t,3),z_da(t,4),obs_chr(t)
    ENDIF
  ENDDO
  DO t=nt_asm+1,nt_asm+nt_prd
    IF(MOD(t,output_interval)==0) THEN
      WRITE(*,'(F7.2,7F10.3)') dt*t,z_t(t),z_s(t,2),z_s(t,3),z_s(t,4),&
        &z_da(t,2),z_da(t,3),z_da(t,4)
    ENDIF
  ENDDO

  WRITE(*,*)
  WRITE(*,*) "---------------------"
  WRITE(*,*) "        RMSE         "
  WRITE(*,*) " [Time] [No DA] [DA] "

  DO t=0,nt_asm+nt_prd
    SQ_s(t,1)=SQ_s(t-1,1)+(x_s(t,3)-x_t(t))**2
    RMSE_s(t,1)=SQRT(SQ_s(t,1)/(t+1))
    SQ_da(t,1)=SQ_da(t-1,1)+(x_da(t,3)-x_t(t))**2
    RMSE_da(t,1)=SQRT(SQ_da(t,1)/(t+1))
    IF (MOD(t,output_interval)==0) THEN
      WRITE(*,'(F7.2,3F10.3)') dt*t, RMSE_s(t,1),RMSE_da(t,1)
    ENDIF
  ENDDO

  STOP

  CONTAINS

  SUBROUTINE tinteg_L63_true(xin,yin,zin,xout,yout,zout)
    IMPLICIT NONE

    REAL(8),PARAMETER :: dt=1.0d-3
    REAL(8),PARAMETER :: parm1=10.0d0,parm2=28.0d0,parm3=8.0d0/3.0d0
    INTEGER,PARAMETER :: n=3
    REAL(8),INTENT(IN) :: xin,yin,zin
    REAL(8),INTENT(OUT) ::xout,yout,zout
    REAL(8) :: M_t(n,n)
    REAL(8) :: X(n)
    REAL(8) :: X3(n)

    X(1)=xin;X(2)=yin;X(3)=zin
    M_t(1,1)=1.0d0-dt*parm1;M_t(1,2)=dt*parm1;M_t(1,3)=0.0d0
    M_t(2,1)=dt*(parm2-X(3));M_t(2,2)=1.0d0-dt;M_t(2,3)=-dt*X(1)
    M_t(3,1)=dt*X(2);M_t(3,2)=dt*X(1);M_t(3,3)=1.0d0-dt*parm3
    X3=MATMUL(M_t,X)
    xout=X3(1);yout=X3(2);zout=X3(3)

    RETURN
  END SUBROUTINE tinteg_L63_true

  SUBROUTINE tinteg_L63_sim(xin,yin,zin,xout,yout,zout,parm1,parm2,parm3)
    IMPLICIT NONE

    REAL(8),PARAMETER :: dt=1.0d-3
    REAL(8),INTENT(IN) :: parm1,parm2,parm3
    INTEGER,PARAMETER :: n=3
    REAL(8),INTENT(IN) :: xin,yin,zin
    REAL(8),INTENT(OUT) ::xout,yout,zout
    REAL(8) :: M_t(n,n)
    REAL(8) :: X(n)
    REAL(8) :: X3(n)

    X(1)=xin;X(2)=yin;X(3)=zin
    M_t(1,1)=1.0d0-dt*parm1;M_t(1,2)=dt*parm1;M_t(1,3)=0.0d0
    M_t(2,1)=dt*(parm2-X(3));M_t(2,2)=1.0d0-dt;M_t(2,3)=-dt*X(1)
    M_t(3,1)=dt*X(2);M_t(3,2)=dt*X(1);M_t(3,3)=1.0d0-dt*parm3
    X3=MATMUL(M_t,X)
    xout=X3(1);yout=X3(2);zout=X3(3)

    RETURN
  END SUBROUTINE tinteg_L63_sim

END PROGRAM Lorenz63_multiSeSRF_DBF
