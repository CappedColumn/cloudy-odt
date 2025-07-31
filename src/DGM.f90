module DGM
    use globals
    use microphysics, only: saturation_vapor_pressure
    implicit none

    ! From DGM-array.f90
    ! Verify where/when these are set and overwritten
    integer(i4), parameter :: nmax = 8
    integer(i4), parameter :: nvar = 8
    
    ! Mani sets these using "set_odeint_dgm_params()"
    integer(i4)            :: dmaxa, ndmax = 1!, err_dgm
    real(dp)               :: grid_scale, solute_mass
    
    ! From DGM-const.f90
    real(dp)               :: Lcond_temp ! Temp. dependent latent heat
    real(dp)               :: Ktemp, D   ! Temp. dependent diffusivities
    real(dp)               :: cpm     ! Specific heat of...

    ! Error tolerance
    real(dp), parameter :: hmin = 0.0

    public :: integrate_ODE, set_aerosol_properties
    private

contains

subroutine set_aerosol_properties(dmax, m0_aerosol, gscale)
    ! Set aerosol properties
    integer(i4), intent(in) :: dmax
    real(dp), intent(in) :: m0_aerosol, gscale

    ! Set the global variables
    dmaxa = dmax
    solute_mass = m0_aerosol
    grid_scale = gscale

end subroutine set_aerosol_properties

! Below is the logic of the ancients, pre-dating most modern forms
! of written human communication. Somewhere, in a creeky basement next
! to a boiler, sits a dusty bookshelf containing the transcription. A
! rosetta stone to illuminate the all-cap intrinsics of yester-years.

! **** Turns on ~THE HISTORY CHANNEL~ **** 
! "The comment was not invented until 1843, by Augustus Chaddington VI, an
! amateur naturalist and tinkerer who's "magic logic machines" were built
! from herring bones, pea-coat buttons and half-combusted pipe tobacco. After
! meeting a lady at his local croquet club, he created the inline comment
! in an attempt to woo her. They married 3 months later."

SUBROUTINE integrate_ODE(ystart,x1,x2,h1)

   INTEGER       :: i, kmax, kmaxx, kount, maxstp
   INTEGER       :: nbad, nok, nstp
   REAL(dp)        :: eps, h1, x1, x2, TINY
   PARAMETER (maxstp=10000,kmaxx=200,TINY=1.e-30)
   REAL(dp)        :: dydx(nmax), ystart(nvar)
   REAL(dp)        :: dxsav, h, hdid, hnext, x, xsav
   REAL(dp)        :: xp(kmaxx), y(nmax), yp(nmax,kmaxx), yscal(nmax)

   eps = 1e-4 ! Error tolerance

   dydx(:) = 0.0

   x     = x1
   h     = SIGN(h1,x2-x1)
   nok   = 0
   nbad  = 0
   kount = 0
   kmax  = 0

   DO i=1,nvar
      y(i) = ystart(i)
   END DO

   IF (kmax .GT. 0) xsav=x-2.*dxsav

   DO nstp=1,maxstp
      CALL fcnkb(x,y,dydx)

      DO i=1,nvar
         yscal(i) = ABS(y(i))+ABS(h*dydx(i))+TINY
      END DO

      IF (kmax .GT. 0) THEN
         IF (ABS(x-xsav) .GT. ABS(dxsav)) THEN
            IF (kount .LT. kmax-1) THEN
               kount     = kount+1
               xp(kount) = x
               DO i=1,nvar
                  yp(i,kount) = y(i)
               END DO
               xsav = x
            END IF
         END IF
      END IF

      IF ((x+h-x2)*(x+h-x1) .GT. 0.0) h=x2-x

      CALL rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext)

      IF (hdid .EQ. h) THEN
         nok  = nok+1
      ELSE
         nbad = nbad+1
      END IF

      IF ((x-x2)*(x2-x1) .GE. 0.0) THEN
         DO i=1,nvar
            ystart(i)=y(i)
         END DO

         IF (kmax .NE. 0) THEN
            kount=kount+1
            xp(kount)=x
            DO i=1,nvar
               yp(i,kount)=y(i)
            END DO
         END IF
         RETURN
      END IF

      IF (ABS(hnext) .LT. hmin) THEN
         WRITE(7,*) 'stepsize smaller than minimum in odeint'
         STOP
      END IF

      h = hnext
   END DO
   WRITE(7,*) 'too many steps in odeint'
   RETURN

   END


subroutine fcnkb(ltime, drop_radius, drdt)
    !This function provides derivates for y variables in drop_radius at ltime
    !input:     ltime               -> time
    !input:     drop_radius     -> y variables
    !output:    drdt            -> time derivative of all y variables
    real(dp), intent(in)    :: ltime, drop_radius(8)
    real(dp), intent(out) :: drdt(8)

    real(dp)  :: ck, cr, denom
    real(dp)  :: e, es
    real(dp)  :: falpha, fbeta, rho, rhol

    ! -- molecular weight of the solute (NaCl, (NH4)2SO4, ...)
    real(dp)  :: Ms

    real(dp)  :: radius, qv, temp, s, press, ver_vel, height, ql, nions

    real(dp)  :: lalpha, lbeta
    real(dp), parameter :: alph = 1
    real(dp), parameter :: beta = 0.04

    real(dp)  :: c7
    real(dp), parameter :: c7am   = 0.4363021
    real(dp), parameter :: c7nacl = 0.5381062
    real(dp), parameter :: sigma  = 7.392730e-2

    nions   = 0.0d0

    radius  = drop_radius(1)
    qv      = drop_radius(2)
    temp    = drop_radius(3)
    s       = drop_radius(4)
    press   = drop_radius(5)
    ver_vel = drop_radius(6)
    height  = drop_radius(7)
    ql      = drop_radius(8)

    !Mani: constanst used from module const
    cpm     = cp*((1.0+cp_wv/cp*qv)/(1.0+qv))

    ! -- initialize diffusivities D, Ktemp and Lcond_temp

    ! -- Equation 2 from Bolton, MWR (1980)
    Lcond_temp      = (2.501-0.00237*(temp-Tice))*1.e6

    ! -- Table 7.1 (p103): linear interpolation between -10 and 30 deg C
    ! -- Rogers and Yau (1989)
    Ktemp   = 7.7e-5*(temp-Tice)+0.02399
    D       = 1.57e-7*(temp-Tice)+2.211e-5
    D       = D*1.e5/press

    es      = saturation_vapor_pressure(temp)

    ! -- calculate lalpha and lbeta
    !Mani: I found that sqrt of negative temp produces floating point error, so I am catching it before occurs
    !if (temp < 0) then
    if (temp < 273 .or. qv < 0 .or. temp > 320) then
      write(*,*) ""
      write(*,'(3(A,E16.8))') "fcnb error: T: ", temp, " qv: ", qv,  " radius: " , radius
      write(*,*) "fcnb error: Ktemp: ", Ktemp
      write(*,*) "fcnb error: 2.0*pi*Md*R_uni : ", 2.0*pi*Ma*1e-3*R_univ
      flush(6)
      stop
      !err_dgm   = -1
      !return
      ! write(*,*) "fcnb error: SQRT(2.0*pi*Md*R_uni*temp): ", SQRT(2.0*pi*Md*R_uni*temp)
      ! flush(6)
    end if

    !Mani: if the temperature oscillates to a negative value, sqrt() will return NaN and throws a floating point 
    !exception. if compiled with Floating point flags on, the exception will terminate the program. 
    !This is difficult to debug. therefore compile without floating point flags, and handle it in the program.
    !So you know, where the issue occurs.

    lalpha  = Ktemp * SQRT(2.0*pi*Ma*R_univ*temp)/(alph*press*(cv+R_univ/2.0))
    lbeta   = SQRT(2.0*pi*Mw/(Rv*temp))*D/beta

    !Mani:dmaxa is used from the module array
    !so set the value before calling odeint
    if (dmaxa == 1) then     
      !sodium cholride NaCl
       Ms = 58.4428e-3 ! kg/m3
       c7 = c7nacl
       nions = 2.0d0
    else if(dmaxa == 2) then
      !Ammonium sulphate
       Ms = 132.1395e-3
       !density of Ammonium bisulphate and sulphate are similar so using the same density factor
       c7 = c7am
       nions = 3.0d0
    else if(dmaxa == 3) then
      !Ammonium bisulphate
        Ms = 115.11e03
        !density of Ammonium bisulphate and sulphate are similar so using the same density factor
        c7 = c7am
        nions = 2.0d0
    else  
      write(*,*) "fcnkb: Error: aerosol type unknown"
      stop
    end if


    falpha  = radius/(radius+lalpha)
    fbeta   = radius/(radius+lbeta)
    rhol    = (radius**3*pi_43*rho_l+solute_mass*c7)/(radius**3*pi_43)

    ! -- calculate drdt
    ck      = (2*sigma/(Rv*temp*rhol*radius))
    cr      = nions*(Mw/Ms)*solute_mass/(pi_43*radius**3*rhol-solute_mass)
    denom   = rhol*(Rv*temp/(fbeta*D*es)+Lcond_temp**2/(falpha*Ktemp*Rv*temp**2))
    drdt(1) = 1.0/radius*(s-ck+cr)/denom
    !write(*,*) "drdt(1): ", drdt(1)
    ! -- calculate change of water vapor
    drdt(2) = -pi_4*grid_scale*radius**2*drdt(1)*rho_l

    ! Verify that water vapor does not go negative
    if ( (drdt(2) < 0.0) .and. (abs(drdt(2)) > qv) ) then
      drdt(2) = - qv
      drdt(1) = - drdt(2)/(pi_4*grid_scale*radius**2*rho_l)
    end if
    
    ! -- calculate change of temperature
    drdt(3) = -Lcond_temp/cpm*drdt(2)   !-ver_vel*g/cpm

    !write(*,*) "drdt(3): ", drdt(3)
    ! -- calculate change of supersaturation
    !e       = qv*press/(eps_err+qv)
    !rho     = (press-e)/(Rd*temp)*(1.0+ql)+e/(Rv*temp)
    !drdt(4) = (1.0+s)*(Rd/(qv*(Rv*qv+Rd))*drdt(2)-g*rho*ver_vel/press-Lcond_temp/(Rv*temp**2)*drdt(3))
    ! -- calculate change of pressure
    !drdt(5) = -1.0*rho*g*ver_vel
    ! -- calculate change of vertical velocity
    !drdt(6) = 0.0
    ! -- calculate change of height
    !drdt(7) = ver_vel
    ! -- calculate change of liquid water
    drdt(8) = -1.0*drdt(2)

end subroutine fcnkb !of costa mesa

SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext)
   INTEGER       :: i, n
   REAL(dp)        :: eps, errcon, errmax, h, hdid, hnext, htry, pgrow, pshrnk, safety, x, xnew
   REAL(dp)        :: dydx(n), y(n), yerr(nmax), yscal(n), ytemp(nmax)
   PARAMETER (safety=0.9,pgrow=-0.2,pshrnk=-0.25,errcon=1.89e-4)

   h      = htry

1       CALL rkck(y,dydx,n,x,h,ytemp,yerr)

   errmax = 0.0
   DO i=1,n
      errmax = MAX(errmax,ABS(yerr(i)/yscal(i)))
   END DO

   errmax = errmax/eps

   IF (errmax .GT. 1.0) THEN
      h = safety*h*(errmax**pshrnk)
      IF (h .LT. 0.1*h) THEN
         h = 0.1*h
      END IF
      xnew = x+h
      IF (xnew .EQ. x) THEN
         WRITE(7,*) 'stepsize underflow in rkqs'
         STOP
      END IF
      GOTO 1
   ELSE
      IF(errmax .GT. errcon) THEN
         hnext = safety*h*(errmax**pgrow)
      ELSE
         hnext = 5.0*h
      END IF
      hdid = h
      x    = x+h
      DO i=1,n
         y(i) = ytemp(i)
    IF ((i .LE. ndmax) .AND. (y(i) .LE. 0.0)) y(i) = 1.d-8
      END DO
      RETURN
   END IF

   END

! -- (C) Copr. 1986-92 Numerical Recipes Software 5"#@-130Rk#3#)KB.

SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr)

   INTEGER       :: i, n
   REAL(dp)        :: dydx(n), h, x, y(n), yerr(n), yout(n)

   REAL(dp)        :: ak2(nmax), ak3(nmax), ak4(nmax), ak5(nmax), ak6(nmax),       &
                 &  ytemp(nmax), A2, A3, A4, A5, A6, B21, B31, B32, B41, B42, B43, B51,  &
                 &  B52, B53, B54, B61, B62, B63, B64, B65, C1, C3, C4, C6, DC1, DC3, DC4, &
                 &  DC5, DC6
   PARAMETER ( A2=0.2,A3=0.3,A4=0.6,A5=1.0,A6=0.875,B21=0.2,B31=3.0/40.0,&
             & B32=9.0/40.0,B41=0.3,B42=-0.9,B43=1.2,B51=-11.0/54.0,B52=2.5,&
             & B53=-70.0/27.0,B54=35.0/27.0,B61=1631.0/55296.0,B62=175.0/512.0,&
             & B63=575.0/13824.0,B64=44275.0/110592.0,B65=253.0/4096.0,C1=37.0/378.0,&
             & C3=250.0/621.0,C4=125.0/594.0,C6=512.0/1771.0,DC1=C1-2825.0/27648.0,&
             & DC3=C3-18575.0/48384.0,DC4=C4-13525.0/55296.0,DC5=-277.0/14336.0,&
             & DC6=C6-0.25 )

   DO i=1,n
      ytemp(i) = y(i)+B21*h*dydx(i)
   END DO

   CALL fcnkb(x+A2*h,ytemp,ak2)

   DO i=1,n
      ytemp(i) = y(i)+h*(B31*dydx(i)+B32*ak2(i))
   END DO

   CALL fcnkb(x+A3*h,ytemp,ak3)

   DO i=1,n
      ytemp(i) = y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
   END DO

   CALL fcnkb(x+A4*h,ytemp,ak4)

   DO i=1,n
      ytemp(i) = y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
   END DO

   CALL fcnkb(x+A5*h,ytemp,ak5)

   DO i=1,n
      ytemp(i) = y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
   END DO

   CALL fcnkb(x+A6*h,ytemp,ak6)

   DO i=1,n
      yout(i)  = y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
   END DO

   DO i=1,n
      yerr(i)  = h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
   END DO
 RETURN

 END

! -- (C) Copr. 1986-92 Numerical Recipes Software 5"#@-130Rk#3#)KB.


end module DGM