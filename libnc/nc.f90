module nc

    implicit none

    private
    public :: check

contains

    subroutine check(istatus)

        ! catch and throw error according to NF90 library

        use netcdf

        integer, intent(in) :: istatus

        if (istatus /= nf90_noerr) then
            write(*,*) trim(adjustl(nf90_strerror(istatus)))
            error stop 1
        end if

    end subroutine check

end module nc