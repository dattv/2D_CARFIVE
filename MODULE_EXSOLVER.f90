MODULE MODULE_EXSOLVER
    
    use MODULE_PRECISION
    use MODULE_CFD_DATA
    use MODULE_QUADTREE
    use MODULE_GENERICMETHOD
    use MODULE_TIMESOLVER
    
    private
    public  :: second_EX
    
    contains
!==================================================================================================
    subroutine Second_EX(first, last, tree, dt)
    implicit none
    integer(ip), intent(in) :: first, last
    real(rp), intent(inout) :: dt
    type(quadtree), dimension(first:last), intent(inout)    :: tree
    
    call loop_on_quadtree_array(first, last, tree, second_Ex_single)
    
    call compute_dt_array(first, last, tree, dt)
    
    return
    end subroutine Second_EX
!==================================================================================================   
    subroutine second_Ex_single(tree)
    implicit none
    type(quadtree), pointer, intent(inout)   :: tree
 
    return
    end subroutine second_EX_single
!================================================================================================== 

!==================================================================================================    
    
END MODULE MODULE_EXSOLVER    