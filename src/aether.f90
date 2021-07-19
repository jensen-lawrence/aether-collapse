! Performs a simulation of spherically symmetric collapse in Einstein-aether theory
! with an additional scalar matter field

! Adapted from code written by Garfinkle et al. for the papers: 
! "Numerical simulations of gravitational collapse in Einstein-aether theory" (https://arxiv.org/abs/gr-qc/0703093)
! "White holes in Einstein-aether theory" (https://arxiv.org/abs/1608.06970)

! ----------------------------------------------------------------------------------------------------------------------
! Program Execution
! ----------------------------------------------------------------------------------------------------------------------
program main 
    implicit none

    ! Initialize variables and parameters
    character(128) :: savedata
    character(128) :: amp1in, amp2in, r0in, sigmain
    character(128) :: rmaxin, ntimein, timein, tmaxin
    character(128) :: c1in, c2in, c3in, c4in

    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: n = 20000
    integer :: itime, ntime, iter, i
    integer :: ichk

    real(kind=dp), parameter :: noflux = 0.0
    real(kind=dp) :: k(n), ar(n), rad(n)
    real(kind=dp) :: dtk(n), dtar(n), dtrad(n)
    real(kind=dp) :: knew(n), arnew(n), radnew(n)
    real(kind=dp) :: dtknew(n), dtarnew(n), dtradnew(n)
    real(kind=dp) :: psi(n), p(n), dtpsi(n), dtp(n)
    real(kind=dp) :: psinew(n), pnew(n), dtpsinew(n), dtpnew(n)
    real(kind=dp) :: grr(n), dtgrr(n), grrnew(n), dtgrrnew(n)
    real(kind=dp) :: rmax, dr, r(n), dt
    real(kind=dp) :: amp1, amp2, sigma, r0
    real(kind=dp) :: cvec(4), c1, c2, c3, c4
    real(kind=dp) :: c13, c14, c123, ck, v
    real(kind=dp) :: cnstr(n), cnstr2(n)
    real(kind=dp) :: alpha(n)
    real(kind=dp) :: w, lam, s1, s2, s3, wt
    real(kind=dp) :: trap(n), trap2(n)
    real(kind=dp) :: curv(n)
    real(kind=dp) :: mass(n), q(n)
    real(kind=dp) :: time
    real(kind=dp) :: cfl, cflmin, tmax
    real(kind=dp) :: tll(n), tnn(n)
    real(kind=dp) :: asq, trapmin, trapmin2

    ! Read in simulation parameter values
    call get_command_argument(1, savedata)
    call get_command_argument(2, amp1in)
    call get_command_argument(3, amp2in)
    call get_command_argument(4, r0in)
    call get_command_argument(5, sigmain)
    call get_command_argument(6, rmaxin)
    call get_command_argument(7, ntimein)
    call get_command_argument(8, timein)
    call get_command_argument(9, tmaxin)
    call get_command_argument(10, c1in)
    call get_command_argument(11, c2in)
    call get_command_argument(12, c3in)
    call get_command_argument(13, c4in)

    ! Set simulation parameter values
    read(amp1in,*) amp1 
    read(amp2in,*) amp2 
    read(r0in,*) r0 
    read(sigmain,*) sigma 
    read(rmaxin,*) rmax 
    read(ntimein,*) ntime 
    read(timein,*) time 
    read(tmaxin,*) tmax 
    read(c1in,*) c1 
    read(c2in,*) c2 
    read(c3in,*) c3 
    read(c4in,*) c4

    ! Set dr and initial dt
    dr = rmax/float(n - 1)
    dt = 0.1*dr

    ! Set cvec values
    cvec(1) = c1
    cvec(2) = c2
    cvec(3) = c3
    cvec(4) = c4
    c13 = c1 + c3
    c14 = c1 + c4
    c123 = c1 + c2 + c3
    ck = -1.0*c14*(1.0 - c13)/c123
    v = sqrt((c14 - 2.0)/(ck*(2.0 + c13 + 3.0*c2)))

    ! Open files for saving results
    open(unit=14, file=trim(savedata)//'_params.txt', status='unknown')
    open(unit=15, file=trim(savedata)//'arfin.txt', status='unknown')
    open(unit=16, file=trim(savedata)//'kfin.txt', status='unknown')
    open(unit=17, file=trim(savedata)//'cnstrfin.txt', status='unknown')
    open(unit=18, file=trim(savedata)//'radfin.txt', status='unknown')
    open(unit=19, file=trim(savedata)//'alphafin.txt', status='unknown')
    open(unit=21, file=trim(savedata)//'trapfin.txt', status='unknown')
    open(unit=22, file=trim(savedata)//'psifin.txt', status='unknown')
    open(unit=23, file=trim(savedata)//'pfin.txt', status='unknown')
    open(unit=24, file=trim(savedata)//'curvfin.txt', status='unknown')
    open(unit=25, file=trim(savedata)//'massfin.txt', status='unknown')
    open(unit=26, file=trim(savedata)//'grrfin.txt', status='unknown')
    open(unit=27, file=trim(savedata)//'kfin2.txt', status='unknown')
    open(unit=28, file=trim(savedata)//'arfin2.txt', status='unknown')
    open(unit=29, file=trim(savedata)//'qfin2.txt', status='unknown')
    open(unit=30, file=trim(savedata)//'trapfin2.txt', status='unknown')
    open(unit=31, file=trim(savedata)//'cnstrfin2.txt', status='unknown')
    open(unit=32, file=trim(savedata)//'asqfin.txt', status='unknown')
    open(unit=33, file=trim(savedata)//'trapmin.txt', status='unknown')
    open(unit=34, file=trim(savedata)//'trapmin2.txt', status='unknown')
    open(unit=35, file=trim(savedata)//'tnnfin.txt', status='unknown')
    open(unit=36, file=trim(savedata)//'attime.txt', status='unknown')
    open(unit=37, file=trim(savedata)//'ttime.txt', status='unknown')
    open(unit=38, file=trim(savedata)//'Psi.txt', status='unknown')
    open(unit=39, file=trim(savedata)//'P.txt', status='unknown')
    open(unit=40, file=trim(savedata)//'ar.txt', status='unknown')
    open(unit=41, file=trim(savedata)//'K.txt', status='unknown')
    open(unit=42, file=trim(savedata)//'alpha.txt', status='unknown')
    open(unit=43, file=trim(savedata)//'Tll.txt', status='unknown')
    open(unit=44, file=trim(savedata)//'Tnn.txt', status='unknown')

    ! Find initial values
    call init(n, psi, p, ar, k, rad, grr, dr, r, amp1, amp2, sigma, r0, cvec)

    ! Loop through time and evolve system
    do itime = 1, ntime
        ichk = itime/200
        ichk = ichk*200

        ! Write psi, p, ar, k, alpha, tll, tnn to data files
        if (ichk == itime) then
            do i = 1, n/8, 10
                write(38,*) time, r(i), psi(i)
                write(39,*) time, r(i), p(i)
                write(40,*) time, r(i), ar(i)
                write(41,*) time, r(i), k(i)
                write(42,*) time, r(i), alpha(i)
                write(43,*) time, r(i), tll(i)
                write(44,*) time, r(i), tnn(i)
            end do

            write(38,*) ''
            write(39,*) ''
            write(40,*) ''
            write(41,*) ''
            write(42,*) ''
            write(43,*) ''
            write(44,*) ''

            write(*,*) itime, time
        end if

        ! Calculate the time derivatives for current variable values
        call evolve(n, psi, p, ar, k, rad, grr, dtpsi, dtp, dtar, dtk, dtrad, &
            & dtgrr, dr, dt, cvec, cnstr, alpha, trap, trap2, curv, mass, q, cnstr2, &
            & tll,tnn)

        do i = 1, n
            psinew(i) = psi(i)
            pnew(i) = p(i)
            arnew(i) = ar(i)
            knew(i) = k(i)
            radnew(i) = rad(i)
            grrnew(i) = grr(i)
        end do
        
        ! Execute Crank-Nicholson iteration
        do iter = 1, 3

            ! Calculate the time derivatives for variable values at next time step
            call evolve(n, psinew, pnew, arnew, knew, radnew, grrnew, dtpsinew, &
                & dtpnew, dtarnew, dtknew, dtradnew, dtgrrnew, dr, dt, cvec, cnstr, alpha, &
                & trap, trap2, curv, mass, q, cnstr2, tll, tnn)

            ! Find the new variable values by averaging the time derivatives
            do i = 2, n-1
                psinew(i) = psi(i) + dt*0.5*(dtpsi(i) + dtpsinew(i))
                pnew(i) = p(i) + dt*0.5*(dtp(i) + dtpnew(i))
                arnew(i) = ar(i) + dt*0.5*(dtar(i) + dtarnew(i))
                knew(i) = k(i) + dt*0.5*(dtk(i) + dtknew(i))
                radnew(i) = rad(i) + dt*0.5*(dtrad(i) + dtradnew(i))
                grrnew(i) = grr(i) + dt*0.5*(dtgrr(i) + dtgrrnew(i))
            end do

            ! Apply the boundary conditions at r = 0 and r = rmax
            psinew(1) = (4.0*psinew(2) - psinew(3))/3.0
            pnew(1) = (4.0*pnew(2) - pnew(3))/3.0
            arnew(1) = 0.0
            knew(1) = (4.0*knew(2) - knew(3))/3.0
            radnew(1) = 0.0
            grrnew(1) = (4.0*grrnew(2) - grrnew(3))/3.0

            psinew(n) = 0.0
            pnew(n) = 0.0
            arnew(n) = arnew(n-1)
            ! knew(n) = 0.0
            knew(n) = noflux * q(n)
            grrnew(n) = 1.0
            radnew(n) = rad(n) + dt*0.5*(dtrad(n) + dtradnew(n))
            grrnew(n) = grr(n) + dt*0.5*(dtgrr(n) + dtgrrnew(n))
        end do

        ! Update current variable values
        do i = 1, n
            psi(i) = psinew(i)
            p(i) = pnew(i)
            ar(i) = arnew(i)
            k(i) = knew(i)
            rad(i) = radnew(i)
            grr(i) = grrnew(i)
        end do

        time = time + dt

        ! Find minimum trapped surface values and write to data files
        trapmin = trap(1)
        trapmin2 = trap2(1)

        do i = 2, n
            if (trap(i) < trapmin) then
                trapmin = trap(i)
            end if

            if (trap2(i) < trapmin2) then
                trapmin2 = trap2(i)
            end if
        end do

        write(33,*) time, trapmin
        write(34,*) time, trapmin2

        ! Implement the Courant condition
        if (time > tmax) goto 88

        cflmin = 1.0
        do i = 1, n
            cfl = sqrt(grr(i))/alpha(i)
            if (cfl < cflmin) then
                cflmin = cfl
            end if
        end do
        cflmin = (1.0/v) * cflmin
        dt = 0.5*dr*cflmin
    end do

    88  continue

    ! Final time has been reached, write variable values to data files
    do i = 1, n
        write(17,*) r(i), cnstr(i)
        write(31,*) r(i), cnstr2(i)
    end do

    do i = 1, n/4
        write(15,*) r(i), ar(i)
        write(16,*) r(i), k(i)
        write(18,*) r(i), rad(i)
        write(19,*) r(i), alpha(i)
        write(21,*) r(i), trap(i)
        write(30,*) r(i), trap2(i)
        write(22,*) r(i), psi(i)
        write(23,*) r(i), p(i)
        asq=ar(i)*ar(i)/grr(i)
        write(32,*) rad(i), asq
        write(35,*) r(i), tnn(i)
        if((rad(i) > 1.5).and.(curv(i) > 0.0)) then
            write(24,*) log(rad(i)), log(curv(i))
        endif
        write(25,*) rad(i), mass(i)
    end do

    do i = 1138, n
        write(26,*) rad(i), grr(i)
        write(27,*) rad(i), k(i)
        write(28,*) rad(i), ar(i)
        write(29,*) rad(i), q(i)
    end do

    ! Write simulation parameters to data file
    write(14,*) 'n: ', n
    write(14,*) 'amp1: ', amp1
    write(14,*) 'amp2: ', amp2
    write(14,*) 'r0: ', r0
    write(14,*) 'sigma: ', sigma
    write(14,*) 'rmax: ', rmax
    write(14,*) 'tmax: ', tmax
    write(14,*) 'ntime: ', ntime
    write(14,*) 'c1: ', c1
    write(14,*) 'c2: ', c2
    write(14,*) 'c3: ', c3
    write(14,*) 'c4: ', c4
    write(14,*) 'v0: ', v
    write(14,*) 'noflux: ', noflux

end program main

! ----------------------------------------------------------------------------------------------------------------------
! Subroutine to Find the Initial Values
! ----------------------------------------------------------------------------------------------------------------------
subroutine init(n, psi, p, ar, k, rad, grr, dr, r, amp1, amp2, sigma, r0, cvec)
    implicit none

    ! Initialize variables and parameters
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: n
    integer :: i 

    real(kind=dp) :: psi(n), p(n), ar(n), k(n), rad(n), grr(n), dr, r(n), q(n)
    real(kind=dp) :: amp1, amp2, sigma, r0, rsq
    real(kind=dp) :: cvec(4)

    ! Set variable values at t = 0
    do i = 2, n 
        r(i) = dr*float(i - 1)
        grr(i) = 1.0
        psi(i) = amp1 * exp(-1.0 * ((r(i)**2 - r0**2)/sigma**2)**2)
        p(i) = 0.0
        ar(i) = 0.0
        k(i) = 0.0
    end do 

    r(1) = 0.0
    grr(1) = 1.0
    psi(1) = (4.0*psi(2) - psi(3))/3.0
    p(1) = (4.0*p(2) - p(3))/3.0
    ar(1) = 1.0
    k(1) = (4.0*k(2) - k(3))/3.0

    ! Solve the Hamiltonian constraint to find the metric component
    call calcqrad(n, q, rad, k, dr, cvec, psi, p, ar)
    return 
end subroutine init

! ----------------------------------------------------------------------------------------------------------------------
! Subroutine to Calculate the Time Derivatives of the Variables
! ----------------------------------------------------------------------------------------------------------------------
subroutine evolve(n, psi, p, ar, k, rad, grr, dtpsi, dtp, dtar, dtk, dtrad, &
    & dtgrr, dr, dt, cvec, cnstr, alpha, trap, trap2, curv, mass, q, cnstr2, &
    & tll,tnn)
    implicit none

    ! Initialize variables and parameters
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: n
    integer :: i 

    real(kind=dp) :: psi(n),p(n),ar(n),k(n),rad(n),grr(n)
    real(kind=dp) :: dtpsi(n),dtp(n),dtar(n),dtk(n),dtrad(n),dtgrr(n)
    real(kind=dp) :: dr,dt
    real(kind=dp) :: drk,drar,drrad,drrrad
    real(kind=dp) :: drpsi,drp,drgrr,drrpsi,lappsi,diva
    real(kind=dp) :: radp,radm,drpsip,drpsim
    real(kind=dp) :: eps,arko(n),kko(n),psiko(n),pko(n),radko(n),grrko(n)
    real(kind=dp) :: alpha(n),trap(n),trap2(n),q(n),cvec(4)
    real(kind=dp) :: c1,c2,c3,c4,c13,c14
    real(kind=dp) :: term1,term2,term3,term4,term5,term6,term7,term8,term9
    real(kind=dp) :: cnstr(n)
    real(kind=dp) :: cof1,cof2,cof3
    real(kind=dp) :: curv(n),mass(n)
    real(kind=dp) :: cnstr2(n),w(n),intw
    real(kind=dp) :: rho,js,s,mss,tll(n),tnn(n)

    ! Unpack values of c parameters 
    c1 = cvec(1)
    c2 = cvec(2)
    c3 = cvec(3)
    c4 = cvec(4)
    c13 = c1 + c3
    c14 = c1 + c4

    ! Compute the Kreiss-Oliger dissipation terms for extra stability
    eps = 0.5
    do i = 1, n
        arko(i) = 0.0
        kko(i) = 0.0
        psiko(i) = 0.0
        pko(i) = 0.0
        radko(i) = 0.0
        grrko(i) = 0.0
    end do

    do i = 3, n-2
        arko(i) = (-1.0/16.0)*(eps/dt)*(ar(i+2) + ar(i-2) + 6.0*ar(i) - 4.0*(ar(i+1) + ar(i-1)))
        kko(i) = (-1.0/16.0)*(eps/dt)*(k(i+2) + k(i-2) + 6.0*k(i) - 4.0*(k(i+1) + k(i-1)))
        psiko(i) = (-1.0/16.0)*(eps/dt)*(psi(i+2) + psi(i-2) + 6.0*psi(i) - 4.0*(psi(i+1) + psi(i-1)))
        pko(i) = (-1.0/16.0)*(eps/dt)*(p(i+2) + p(i-2) + 6.0*p(i) - 4.0*(p(i+1) + p(i-1)))
        radko(i) = (-1.0/16.0)*(eps/dt)*(rad(i+2) + rad(i-2) + 6.0*rad(i) - 4.0*(rad(i+1) + rad(i-1)))
        grrko(i) = (-1.0/16.0)*(eps/dt)*(grr(i+2) + grr(i-2) + 6.0*grr(i) - 4.0*(grr(i+1) + grr(i-1)))
    end do

    ! Find the lapse, trace-free part of extrinsic curvature, and shift
    call calcalpha(n, alpha, ar, dr)
    call calcq(n, q, rad, k, dr, cvec, psi, p)

    ! Calculate the derivatives of the variables
    do i = 2, n-1
        ! Spatial derivatives
        drk = (k(i+1) - k(i-1))/(2.0*dr)
        drar = (ar(i+1) - ar(i-1))/(2.0*dr)
        radp = rad(i+1)
        radm = rad(i-1)
        diva = (1.5/(dr*(radp*radp + radp*radm + radm*radm)))*(radp*radp*ar(i+1) - radm*radm*ar(i-1))
        drrad = (rad(i+1) - rad(i-1))/(2.0*dr)
        drpsi = (psi(i+1) - psi(i-1))/(2.0*dr)
        drp = (p(i+1) - p(i-1))/(2.0*dr)
        drgrr = (grr(i+1) - grr(i-1))/(2.0*dr)

        drrpsi = (psi(i+1) + psi(i-1) - 2.0*psi(i))/dr**2
        drrrad = (rad(i+1) + rad(i-1) - 2.0*rad(i))/dr**2
        radp = 0.5*(rad(i+1) + rad(i))
        radm = 0.5*(rad(i-1) + rad(i))
        drpsip = (psi(i+1) - psi(i))/dr
        drpsim = (psi(i) - psi(i-1))/dr
        lappsi = (3.0/(dr*(radp*radp + radp*radm + radm*radm)))*(radp*radp*drpsip - radm*radm*drpsim)

        ! Time derivatives
        term1 = alpha(i)*p(i)
        term2 = psiko(i)
        dtpsi(i) = term1 + term2

        term1 = alpha(i)*(p(i)*k(i) + ar(i)*drpsi/grr(i))
        term2 = alpha(i)*lappsi/grr(i)
        term3 = -0.5*alpha(i)*drgrr*drpsi/(grr(i)*grr(i))
        term4 = pko(i)
        dtp(i) = term1 + term2 + term3 + term4
        cof1 = c13/(c14*(1.0 - c13))
        cof2 = -1.0*(c2 + c13)/(c14*(1.0 - c13))

        term1 = alpha(i)*(k(i)/3.0 - 2.0*q(i))*ar(i)
        term2 = alpha(i)*cof1*p(i)*drpsi
        term3 = alpha(i)*cof2*drk
        term4 = arko(i)
        dtar(i) = term1 + term2 + term3 + term4
        cof1 = (c14 - 2.0)/(2.0 + c13 + 3.0*c2)
        cof2 = 2.0/(2.0 + c13 + 3.0*c2)
        cof3 = 3.0*(1.0 - c13)/(2.0 + c13 + 3.0*c2)

        term1 = (alpha(i)/3.0)*k(i)*k(i)
        term2 = alpha(i)*cof1*(diva + ar(i)*ar(i))/grr(i)
        term3 = -0.5*alpha(i)*cof1*drgrr*ar(i)/(grr(i)*grr(i))   
        term4 = alpha(i)*cof2*p(i)*p(i)
        term5 = alpha(i)*cof3*q(i)*q(i)
        term6 = kko(i)
        dtk(i) = term1 + term2 + term3 + term4 + term5 + term6

        term1 = alpha(i)*rad(i)*(q(i)/2.0 - k(i)/3.0)
        term2 = radko(i)
        dtrad(i) = term1 + term2
        ! modify ko term?

        term1 = -2.0*alpha(i)*grr(i)*(q(i) + k(i)/3.0)
        term2 = grrko(i)
        dtgrr(i) = term1 + term2

        ! Calculate the first constraint to validate results
        term1 = drrrad
        term2 = (0.5/rad(i))*(drrad**2 - grr(i))
        term3 = (-0.5/grr(i))*drgrr*drrad
        term4 = c14*rad(i)*(0.5*drar + 0.25*ar(i)**2)
        term5 = c14*ar(i)*drrad
        term6 = -0.25*(rad(i)*c14/grr(i))*drgrr*ar(i)
        term7 = (3.0/8.0)*rad(i)*grr(i)*(1.0 - c13)*(q(i)**2)  
        term8 = (-1.0/12.0)*rad(i)*grr(i)*(2.0 + c13 + 3.0*c2)*(k(i)**2)  
        term9 = 0.25*rad(i)*(grr(i)*(p(i)**2)+drpsi**2)
        cnstr(i) = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9
        w(i) = (2.0*rad(i)*drrad/grr(i))*(term4 + term5 + term6 + term7 + term8 + term9)
        cnstr(i) = cnstr(i)*rad(i)/(10.0 + rad(i))

        ! Check for trapped and anti-trapped surfaces and calculate mass
        trap(i) = rad(i)*(q(i)/2.0 - k(i)/3.0) + drrad/sqrt(grr(i))
        trap2(i) = -1.0*rad(i)*(q(i)/2.0 - k(i)/3.0) + drrad/sqrt(grr(i))
        mass(i) = (rad(i)/2.0)*(1.0 - rad(i))

        ! Calculate the stress-energy tensor components
        term1 = (0.5*c14/grr(i))*(2.0*drar + ar(i)*(4.0*drrad/rad(i) + ar(i)))
        term2 = 0.5*(p(i)*p(i) + drpsi*drpsi/grr(i))
        term3 = (-1.0/6.0)*(c13 + 3.0*c2)*k(i)*k(i)
        term4 = -0.75*c13*q(i)*q(i)
        term5 = -0.5*c14*ar(i)*drgrr/(grr(i)*grr(i))
        rho = term1 + term2 + term3 + term4 + term5
        js = (1.0/(1.0 - c13))*((c13 + c2)*drk - p(i)*drpsi)/(sqrt(grr(i)))

        term1 = 2.0*p(i)*p(i) - 4.5*(c13 + c2)*q(i)*q(i)
        term2 = (c14 + c13 + 2.0*c2)*(drar + ar(i)*(2.0*drrad/rad(i) + ar(i)))/grr(i)
        term3 = -0.5*(c14 + c13 + 2.0*c2)*ar(i)*drgrr/(grr(i)*grr(i))
        s = -1.0*rho + (2.0/(2.0 + c13 + 3.0*c2))*(term1 + term2 + term3)

        term1 = drpsi*drpsi + (c3 - c4)*ar(i)*ar(i)
        term2 = c13*(drar + ar(i)*drrad/rad(i))
        term3 = (c13/rad(i))*(drrrad + (grr(i) - drrad*drrad)/rad(i))
        term4 = (-0.5*drgrr/grr(i))*(c13*ar(i) + drrad/rad(i))
        mss = (2.0*grr(i)/(3.0*(1.0 - c13)))*(term1 + term2 + term3 + term4)
        tll(i) = rho + (s/3.0) + mss - 2.0*js
        tnn(i) = rho + (s/3.0) + mss + 2.0*js
    end do

    ! Set values at r = 0 and r = rmax
    cnstr(1) = cnstr(2)
    cnstr(n-1) = cnstr(n-2)
    cnstr(n) = cnstr(n-1)
    trap(1) = 1.0
    trap(n) = trap(n-1)
    trap2(1) = sqrt(2.0)
    trap2(n) = trap2(n-1)
    mass(1) = 0.0
    mass(n) = mass(n-1)
    dtrad(n) = alpha(n)*rad(n)*(q(n)/2.0 - k(n)/3.0)
    dtgrr(n) = -2.0*alpha(n)*grr(n)*(q(n) + k(n)/3.0)

    ! Calculate the second constraint to validate results
    cnstr2(1) = 0.0
    w(1) = 0.0
    intw = 0.0

    do i = 2, n-2
        drrad = (rad(i+1) - rad(i-1))/(2.0*dr)
        intw = intw + 0.5*dr*(w(i) + w(i-1))
        cnstr2(i) = rad(i)*(1.0 - drrad*drrad/grr(i)) - intw
    end do 

    cnstr2(n-1) = 0.0
    cnstr2(n) = 0.0

    ! Find the squared curvature; blows up at a singularity
    call calccurv(n, curv, ar, k, q, rad, psi, p, dr, cvec)
    return
end subroutine evolve

! ----------------------------------------------------------------------------------------------------------------------
! Subroutine to Calculate the Squared Curvature
! ----------------------------------------------------------------------------------------------------------------------
subroutine calccurv(n, curv, ar, k, q, rad, psi, p, dr, cvec)
    implicit none

    ! Initialize variables and parameters
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: n
    integer :: i 

    real(kind=dp) :: curv(n), ar(n), k(n), q(n), rad(n), psi(n), p(n), dr, cvec(4)
    real(kind=dp) :: c1, c2, c3, c4, c13, c14
    real(kind=dp) :: y, r1, r2, k1, k2, l1, l2, m1, m2
    real(kind=dp) :: term1, term2, term3, term4
    real(kind=dp) :: drk, drar, drrad, drpsi, drrrad 
    real(kind=dp) :: luk

    ! Unpack values of c parameters 
    c1 = cvec(1)
    c2 = cvec(2)
    c3 = cvec(3)
    c4 = cvec(4)
    c13 = c1 + c3
    c14 = c1 + c4

    ! Find the squared curvature
    do i = 2, n-1
        drk = (k(i+1) - k(i-1))/(2.0*dr)
        drar = (ar(i+1) - ar(i-1))/(2.0*dr)
        drrad = (rad(i+1) - rad(i-1))/(2.0*dr)
        drpsi = (psi(i+1) - psi(i-1))/(2.0*dr)
        drrrad = (rad(i+1) + rad(i-1) - 2.0*rad(i))
        y = (p(i)*drpsi - (c13+c2)*drk)/(1.0 - c13)
        r1 = -2.0*drrrad/rad(i)
        r2 = (1.0 - (rad(i)*drrrad + drrad*drrad))/(rad(i)*rad(i))
        k1 = q(i) + (k(i)/3.0)
        k2 = (k(i)/3.0) - (q(i)/2.0)
        l1 = r1 + 2.0*k1*k2
        l2 = r2 + k2*(k1 + k2)

        term1 = (k(i)*k(i))/3.0
        term2 = (c14 - 2.0)*(drar + ar(i)*(ar(i) + 2.0*drrad/rad(i)))
        term3 = 2.0*p(i)*p(i)
        term4 = 3.0*(1.0 - c13)*q(i)*q(i)
        luk = term1 + (term2 + term3 + term4)/(2.0 + c13 + 3.0*c2) 

        term1 = -1.0*l2
        term2 = 0.5*(c1 + c4)*(drar + ar(i)*(ar(i) + 2.0*drrad/rad(i)))  
        term3 = 0.5*(c13 + c2)*(luk - k(i)*k(i))
        term4 = c13*(k(i)*k2 + ar(i)*drrad/rad(i))
        m2 = (term1 + term2 + term3 + term4)/(c13 - 1.0)

        term1 = (drar + ar(i)*(ar(i) + 2.0*drrad/rad(i)))
        term2 = luk - 2.0*m2
        m1 = term1 + term2

        term1 = 3.0*l1*l1 + 4.0*l2*(l2 - l1)
        term2 = 4.0*m1*m1 - 8.0*m2*m2
        term3 = -4.0*y*y
        curv(i) = term1 + term2 + term3
    end do

    curv(1)=curv(2)
    curv(n)=curv(n-1)
    return
end subroutine calccurv

! ----------------------------------------------------------------------------------------------------------------------
! Subroutine to Calculate the Lapse by Integrating the Acceleration of the Aether Field
! ----------------------------------------------------------------------------------------------------------------------
subroutine calcalpha(n, alpha, ar, dr)
    implicit none

    ! Initialize variables and parameters
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: n
    integer :: i 

    real(kind=dp) :: alpha(n),ar(n),dr
    real(kind=dp) :: lnalpha(n)

    ! Calculate the lapse
    lnalpha(n) = 0.0

    do i = n, 2, -1
        lnalpha(i-1) = lnalpha(i) - dr*0.5*(ar(i) + ar(i-1))
    end do

    do i = 1, n
        alpha(i) = exp(lnalpha(i))
    end do
    return
end subroutine calcalpha

! ----------------------------------------------------------------------------------------------------------------------
! Subroutine to Calculate the Trace-Free Part of the Extrinsic Curvature by Integrating the Momentum Constraint
! ----------------------------------------------------------------------------------------------------------------------
subroutine calcq(n, q, rad, k, dr, cvec, psi, p)
    implicit none

    ! Initialize variables and parameters
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: n
    integer :: i

    real(kind=dp) :: q(n), rad(n), k(n), dr, cvec(4), psi(n), p(n)
    real(kind=dp) :: c2, c13
    real(kind=dp) :: radavg, drrad, drk, qint
    real(kind=dp) :: cof, kavg
    real(kind=dp) :: drpsi, pavg
    real(kind=dp) :: term1, term2

    ! Unpack values of c parameters 
    c2 = cvec(2)
    c13 = cvec(1) + cvec(3)
    cof = (2.0 + c13 + 3.0*c2)/(3.0 - 3.0*c13)

    ! Calculate the trace-free part of the extrinsic curvature
    q(1) = 0.0
    do i = 1, n-1
        radavg = 0.5*(rad(i) + rad(i+1))
        kavg = 0.5*(k(i) + k(i+1))
        pavg = 0.5*(p(i) + p(i+1))
        drrad = (rad(i+1) - rad(i))/dr
        drk = (k(i+1) - k(i))/dr
        drpsi = (psi(i+1) - psi(i))/dr
        term1 = (-3.0*q(i)/radavg)*drrad
        term2 = cof*drk - pavg*drpsi/(1.0 - c13) 
        qint = q(i) + 0.5*dr*(term1 + term2)
        term1 = (-3.0*qint/radavg)*drrad
        q(i+1) = q(i) + dr*(term1 + term2)
    end do
    return
end subroutine calcq

! ----------------------------------------------------------------------------------------------------------------------
! Subroutine to Calculate the Trace-Free Part of the Extrinsic Curvature by Integrating the Momentum Constraint
! ----------------------------------------------------------------------------------------------------------------------
subroutine calcqrad(n, q, rad, k, dr, cvec, psi, p, ar)
    implicit none

    ! Initialize variables and parameters
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: n
    integer :: i

    real*8 q(n), rad(n), k(n), dr, cvec(4), psi(n), p(n), ar(n)
    real*8 c2, c13, c14
    real*8 drk, qint, w, wint, radint
    real*8 cof, kavg
    real*8 drpsi, pavg, drar, aravg, s1, s2
    real*8 term1, term2, term3, term4

    ! Unpack values of c parameters
    c2 = cvec(2)
    c13 = cvec(1)+cvec(3)
    c14 = cvec(1)+cvec(4)
    cof = (2.0 + c13 + 3.0*c2)/3.0

    ! Calculate the trace-free part of the extrinsic curvature
    q(1) = 0.0
    rad(1) = 0.0
    q(2) = (0.4/(1.0 - c13))*(cof*(k(2) - k(1)) - p(1)*(psi(2) - psi(1)))
    rad(2) = dr
    w = 1.0 - (dr/12.0)*(6.0*c14*ar(2) + (p(1)*p(1) - cof*k(1)*k(1))*dr)

    do i = 2, n-1
        kavg = 0.5*(k(i) + k(i+1))
        pavg = 0.5*(p(i) + p(i+1))
        aravg = 0.5*(ar(i) + ar(i+1))
        drk = (k(i+1) - k(i))/dr
        drpsi = (psi(i+1) - psi(i))/dr
        drar = (ar(i+1) - ar(i))/dr
        s1 = (cof*drk - pavg*drpsi)/(1. - c13)

        term1 = c14*(2.0*drar + aravg*aravg)
        term2 = pavg*pavg + drpsi*drpsi
        term3 = -1.0*cof*kavg*kavg
        s2 = 0.25*(term1 + term2 + term3)
        radint = rad(i) + 0.5*dr*w
        qint = q(i) + 0.5*dr*(s1 - 3.0*q(i)*w/rad(i))

        term1 = (w*w - 1.0)/(2.0*rad(i))
        term2 = c14*aravg*w
        term3 = rad(i)*s2
        term4 = (3.0/8.0)*(1.0 - c13)*rad(i)*q(i)*q(i)
        wint = w - 0.5*dr*(term1 + term2 + term3 + term4)
        rad(i+1) = rad(i) + dr*wint
        q(i+1) = q(i) + dr*(s1-3.0*qint*wint/radint)

        term1 = (wint*wint - 1.0)/(2.0*radint)
        term2 = c14*aravg*wint
        term3 = radint*s2
        term4 = (3.0/8.0)*(1.0 - c13)*radint*qint*qint
        w = w - dr*(term1 + term2 + term3 + term4)
    end do
    return
end subroutine calcqrad

! ----------------------------------------------------------------------------------------------------------------------