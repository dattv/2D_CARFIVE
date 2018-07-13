MODULE MODULE_EXSOLVER
    
    use MODULE_PRECISION
    use MODULE_CFD_DATA
    use MODULE_QUADTREE
    
    contains
!==================================================================================================
    subroutine Second_EX(first, last, tree, dt)
    implicit none
    integer(ip), intent(in) :: first, last
    real(rp), intent(out)   :: dt
    type(quadtree), dimension(first:last), intent(inout)    :: tree
    
    return
    end subroutine Second_EX
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
    
END MODULE MODULE_EXSOLVER    