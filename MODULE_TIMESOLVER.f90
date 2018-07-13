MODULE MODULE_TIMESOLVER
    
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    use MODULE_QUADTREE
    
    contains
!=================================================================================================    
    subroutine compute_dt_array(first, last, tree, dt)
    implicit none
    integer(ip), intent(in)     :: first, last
    type(quadtree), dimensioN(first:last), intent(in)   :: tree
    real(rp), intent(inout)     :: dt
    integer(ip)                 :: i
    
    do i = first, last
        call compute_dt_single(tree(i), dt)
    end do
    
    return
    end subroutine compute_dt_array
!==================================================================================================  
    recursive subroutine compute_dt_single(tree, dt)
    implicit none
    type(quadtree), intent(in)  :: tree
    real(rp), intent(inout)     :: dt
    
    if (.not. tree%is_leaf) then 
        call compute_dt_single(tree%north_west, dt)
        call compute_dt_single(tree%north_east, dt)
        call compute_dt_single(tree%south_west, dt)
        call compute_dt_single(tree%south_east, dt)
    else
        dt = min(dt, tree%data%dt)
    end if
    return
    end subroutine compute_dt_single
END MODULE MODULE_TIMESOLVER    