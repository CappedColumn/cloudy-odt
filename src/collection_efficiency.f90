module collection_efficiency
    use globals, only: dp, pi
    implicit none
    private

    character(32), public :: coalescence_kernel = 'long'

    integer, parameter :: KERN_LONG = 1, KERN_HALL = 2, KERN_UNITY = 3
    integer :: ikernel = KERN_LONG

    public :: set_kernel_selector, collection_efficiency_E

    ! Hall (1980) table, corrected by SKK (2009, 2017, 2018)
    ! Collector radii in ascending order (μm)
    integer, parameter :: NR = 12, NP = 21

    real(dp), parameter :: r_hall(NR) = &
        [5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, &
         100.0, 150.0, 200.0, 300.0]

    real(dp), parameter :: p_hall(NP) = &
        [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, &
         0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, &
         0.80, 0.85, 0.90, 0.95, 1.00]

    ! E_hall(ir, ip) — rows = collector radius (ascending), cols = ratio
    ! Row 1: r0=5 μm (set equal to 10 μm per SKK 2018)
    ! Row 12: r0=300 μm
    real(dp), parameter :: E_hall(NR, NP) = reshape([ &
        ! p=0.00
        1e-4,  1e-4,  1e-4,  1e-4,  0.001, 0.005, 0.05,  0.20, &
        0.50,  0.77,  0.87,  0.97, &
        ! p=0.05
        1e-4,  1e-4,  1e-4,  1e-4,  0.001, 0.005, 0.05,  0.20, &
        0.50,  0.77,  0.87,  0.97, &
        ! p=0.10
        1e-4,  1e-4,  1e-4,  0.002, 0.07,  0.40,  0.43,  0.58, &
        0.79,  0.93,  0.96,  1.00, &
        ! p=0.15
        1e-4,  1e-4,  0.005, 0.02,  0.28,  0.60,  0.64,  0.75, &
        0.91,  0.97,  0.98,  1.00, &
        ! p=0.20
        0.014, 0.014, 0.016, 0.04,  0.50,  0.70,  0.77,  0.84, &
        0.95,  0.97,  1.00,  1.00, &
        ! p=0.25
        0.017, 0.017, 0.022, 0.085, 0.62,  0.78,  0.84,  0.88, &
        0.95,  1.00,  1.00,  1.00, &
        ! p=0.30
        0.019, 0.019, 0.030, 0.17,  0.68,  0.83,  0.87,  0.90, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.35
        0.022, 0.022, 0.043, 0.27,  0.74,  0.86,  0.89,  0.92, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.40
        0.027, 0.027, 0.052, 0.40,  0.78,  0.88,  0.90,  0.94, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.45
        0.030, 0.030, 0.064, 0.50,  0.80,  0.90,  0.91,  0.95, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.50
        0.033, 0.033, 0.072, 0.55,  0.80,  0.90,  0.91,  0.95, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.55
        0.035, 0.035, 0.079, 0.58,  0.80,  0.90,  0.91,  0.95, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.60
        0.037, 0.037, 0.082, 0.59,  0.78,  0.90,  0.91,  0.95, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.65
        0.038, 0.038, 0.080, 0.58,  0.77,  0.89,  0.91,  0.95, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.70
        0.038, 0.038, 0.076, 0.54,  0.76,  0.88,  0.92,  0.95, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.75
        0.037, 0.037, 0.067, 0.51,  0.77,  0.88,  0.93,  0.97, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.80
        0.036, 0.036, 0.057, 0.49,  0.77,  0.89,  0.95,  1.00, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.85
        0.035, 0.035, 0.048, 0.47,  0.78,  0.92,  1.00,  1.02, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.90
        0.032, 0.032, 0.040, 0.45,  0.79,  1.01,  1.03,  1.04, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=0.95
        0.029, 0.029, 0.033, 0.47,  0.95,  1.30,  1.70,  2.30, &
        1.00,  1.00,  1.00,  1.00, &
        ! p=1.00
        0.027, 0.027, 0.027, 0.52,  1.40,  2.30,  3.00,  4.00, &
        1.00,  1.00,  1.00,  1.00  &
    ], [NR, NP])

contains

    subroutine set_kernel_selector()
        select case (trim(coalescence_kernel))
            case ('long');  ikernel = KERN_LONG
            case ('hall');  ikernel = KERN_HALL
            case ('unity'); ikernel = KERN_UNITY
            case default;   ikernel = KERN_LONG
        end select
    end subroutine set_kernel_selector


    pure function collection_efficiency_E(r_large, r_small) result(E)
        real(dp), intent(in) :: r_large, r_small
        real(dp) :: E

        select case (ikernel)
            case (KERN_LONG);  E = long_kernel(r_large, r_small)
            case (KERN_HALL);  E = hall_kernel(r_large, r_small)
            case (KERN_UNITY); E = 1.0
            case default;      E = long_kernel(r_large, r_small)
        end select
    end function collection_efficiency_E


    pure function long_kernel(r_collector, r_collectee) result(E)
        real(dp), intent(in) :: r_collector, r_collectee
        real(dp) :: E, R_um, r_um, p

        R_um = r_collector * 1.0e6
        r_um = r_collectee * 1.0e6
        p = r_um / max(R_um, 1.0e-3)

        if (R_um < 50.0) then
            E = 4.5e-4 * R_um * R_um * (1.0 - 3.0 / (max(R_um, 1.0e-3) * max(p, 0.01)))
        else
            E = 1.0 - 0.5 * exp(-0.003 * R_um * p)
        end if

        E = max(0.0, min(1.0, E))
    end function long_kernel


    pure function hall_kernel(r_collector, r_collectee) result(E)
        real(dp), intent(in) :: r_collector, r_collectee
        real(dp) :: E, R_um, p

        R_um = r_collector * 1.0e6
        p = r_collectee / max(r_collector, 1.0e-30)
        p = max(0.0, min(1.0, p))

        E = bilinear_interp(R_um, p)
        E = max(0.0, E)
    end function hall_kernel


    pure function bilinear_interp(x, y) result(z)
        real(dp), intent(in) :: x, y
        real(dp) :: z, xc, yc, t, u
        integer :: ix, iy

        xc = max(r_hall(1), min(x, r_hall(NR)))
        yc = max(p_hall(1), min(y, p_hall(NP)))

        ix = NR - 1
        do ix = 1, NR - 1
            if (r_hall(ix+1) >= xc) exit
        end do

        iy = NP - 1
        do iy = 1, NP - 1
            if (p_hall(iy+1) >= yc) exit
        end do

        t = (xc - r_hall(ix)) / max(r_hall(ix+1) - r_hall(ix), 1.0e-30)
        u = (yc - p_hall(iy)) / max(p_hall(iy+1) - p_hall(iy), 1.0e-30)

        z = (1.0 - t) * (1.0 - u) * E_hall(ix, iy)   &
          + t         * (1.0 - u) * E_hall(ix+1, iy)   &
          + (1.0 - t) * u         * E_hall(ix, iy+1)   &
          + t         * u         * E_hall(ix+1, iy+1)
    end function bilinear_interp

end module collection_efficiency
