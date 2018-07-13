MODULE MODULE_RESIDUAL
    
    use MODULE_PRECISION
    use MODULE_QUADTREE
    use MODULE_CONSTANTS
    
    contains
!================================================================================================== 
    subroutine residual_norm(first, last, tree, res_norm)
    implicit none
    integer(ip), intent(in)                                 :: first, last
    type(quadtree), dimension(first:last), intent(in)       :: tree
    real(rp), dimension(tree(first)%data%NQ), intent(out)   :: res_norm
    
    return
    end subroutine residual_norm
!==================================================================================================    
END MODULE MODULE_RESIDUAL    