module interp

    implicit none

    private
    public :: interp1
    
contains
    
    subroutine interp1(yinterp, xv, yv, x, ydefault)

        ! linear interpolation subroutine
    
        double precision, dimension(:), intent(out) :: yinterp
        double precision, dimension(:), intent(in) :: xv, yv, x
        double precision, intent(in) :: ydefault
        integer :: nrowsinterp, nrowsdata, i, j
        logical :: is_ascending
    
        nrowsinterp = size(x)
        nrowsdata = size(xv)
    
        is_ascending = xv(1) <= xv(nrowsdata)

        do i=1,nrowsinterp
            if (is_ascending) then
                if ((x(i) < xv(1)) .or. (x(i) > xv(nrowsdata))) then
                    yinterp(i) = ydefault
                else
                    do j=2,nrowsdata
                        if (x(i) <= xv(j)) then
                            yinterp(i) = (x(i)-xv(j-1)) / (xv(j)-xv(j-1)) * (yv(j)-yv(j-1)) + yv(j-1)
                            exit
                        end if
                    enddo
                end if
            else
                if ((x(i) > xv(1)) .or. (x(i) < xv(nrowsdata))) then
                    yinterp(i) = ydefault
                else
                    do j=2,nrowsdata
                        if (x(i) >= xv(j)) then
                            yinterp(i) = (x(i)-xv(j-1)) / (xv(j)-xv(j-1)) * (yv(j)-yv(j-1)) + yv(j-1)
                            exit
                        end if
                    enddo
                end if
            end if
        enddo
    
    
    end subroutine interp1
    
    
end module interp
    