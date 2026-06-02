module entrainment
    use globals
    use microphysics, only: virtual_temp, update_supersat
    implicit none

    private
    public :: initialize_entrainment, apply_entrainment
    public :: ent_rate, n_blob, psigma, random_entrainment

    ! --- ENTRAINMENT namelist variables ---
    real(dp) :: ent_rate = 2.0
    integer(i4) :: n_blob = 1
    real(dp) :: psigma = 0.1
    logical  :: random_entrainment = .true.

    ! --- Entrainment timing ---
    real(dp) :: t_next_entrain = 0.0

contains

    subroutine initialize_entrainment(vel)
        real(dp), intent(in) :: vel
        integer :: nml_unit, ierr
        character(256) :: nml_line, io_emsg

        namelist /ENTRAINMENT/ ent_rate, n_blob, psigma, random_entrainment

        write(*,*) 'Reading ENTRAINMENT namelist values...'
        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, &
             action='read', status='old')
        if (ierr /= 0) then
            write(0,*) io_emsg; stop 1
        end if
        read(nml=ENTRAINMENT, unit=nml_unit, iostat=ierr)
        if (ierr /= 0) then
            backspace(nml_unit)
            read(nml_unit,'(a)') nml_line
            write(0,'(a)') 'Invalid ENTRAINMENT namelist parameter: '//trim(nml_line)
            stop 1
        end if
        close(nml_unit)

        if (psigma <= 0.0 .or. psigma >= 1.0) then
            write(0,*) 'Error: psigma must be in (0, 1), got: ', psigma
            stop 1
        end if
        if (psigma * n_blob >= 1.0) then
            write(0,*) 'Error: psigma * n_blob must be < 1'
            stop 1
        end if
        if (ent_rate <= 0.0) then
            write(0,*) 'Error: ent_rate must be > 0, got: ', ent_rate
            stop 1
        end if

        t_next_entrain = compute_dt_entm(vel)

        write(*,*) 'entrainment:       ON'
        write(*,*) '  ent_rate:        ', ent_rate
        write(*,*) '  n_blob:          ', n_blob
        write(*,*) '  psigma:          ', psigma
        write(*,*) '  random_entrain:  ', random_entrainment
        write(*,*) '  first dt_entm:   ', t_next_entrain, ' s'

    end subroutine initialize_entrainment


    subroutine apply_entrainment(T_env, qv_env, vel)
        real(dp), intent(in) :: T_env, qv_env, vel
        integer, parameter :: max_blobs = 10
        integer :: blob_start(max_blobs), blob_end(max_blobs), n_final
        integer :: i, k

        if (abs(vel) < 1.0e-30) return
        if (time < t_next_entrain) return

        call place_blobs(blob_start, blob_end, n_final)

        do i = 1, n_final
            do k = blob_start(i), blob_end(i)
                T(k) = T_env
                WV(k) = qv_env
            end do
        end do

        do k = 1, N
            Tv(k) = virtual_temp(T(k), WV(k))
        end do
        call update_supersat(T, WV, SS, pres)

        t_next_entrain = time + compute_dt_entm(vel)

    end subroutine apply_entrainment


    function compute_dt_entm(vel) result(dt_entm)
        real(dp), intent(in) :: vel
        real(dp) :: dt_entm, u

        dt_entm = (real(n_blob, dp) / ent_rate) * (psigma / (1.0 - psigma)) &
                  / abs(vel)

        if (random_entrainment) then
            call random_number(u)
            dt_entm = dt_entm * (-log(1.0 - u))
        end if
    end function compute_dt_entm


    subroutine place_blobs(starts, ends, n_final)
        integer, intent(out) :: starts(:), ends(:), n_final
        integer :: xn, blob_size, s, i
        real(dp) :: u

        blob_size = int(psigma * N)
        xn = N - blob_size * n_blob

        call random_number(u)
        s = int(u * xn) + 1

        n_final = 0
        do i = 1, n_blob
            if (i > 1) s = ends(n_final) + int(real(xn, dp) / n_blob) + 1

            if (s <= N .and. s + blob_size - 1 > N) then
                n_final = n_final + 1
                starts(n_final) = s
                ends(n_final) = N
                n_final = n_final + 1
                starts(n_final) = 1
                ends(n_final) = blob_size - (N - s + 1)
            else if (s > N) then
                n_final = n_final + 1
                starts(n_final) = mod(s - 1, N) + 1
                ends(n_final) = mod(s - 1 + blob_size - 1, N) + 1
                if (ends(n_final) < starts(n_final)) then
                    ends(n_final) = N
                    n_final = n_final + 1
                    starts(n_final) = 1
                    ends(n_final) = blob_size - (N - starts(n_final - 1) + 1)
                end if
            else
                n_final = n_final + 1
                starts(n_final) = s
                ends(n_final) = s + blob_size - 1
            end if
        end do
    end subroutine place_blobs

end module entrainment
