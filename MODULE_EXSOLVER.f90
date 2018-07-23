!=================================================================================================
!> CARTESIAN ADAPTIVE FIVE EQUATION MODEL
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/2D_CARFIVE
!================================================================================================= 
MODULE MODULE_EXSOLVER
    
    use MODULE_PRECISION
    use MODULE_CFD_DATA
    use MODULE_QUADTREE
    use MODULE_GENERICMETHOD
    use MODULE_TIMESOLVER
    use MODULE_MUSCL
    use MODULE_CFDMAINDATA
    use MODULE_RHS
    
    private
    public  :: second_EX
    
    contains
!==================================================================================================
    subroutine Second_EX(first, last, tree, dt)
    implicit none
    integer(ip), intent(in) :: first, last
    real(rp), intent(inout) :: dt
    type(quadtree), dimension(first:last), intent(inout)    :: tree
    
    
    ! ===> muscl reconstruction
    call muscle(int(one  ),  I_limiter_type, first, last, tree)
    call muscle(int(two  ),  I_limiter_type, first, last, tree)
    call muscle(int(three),  I_limiter_type, first, last, tree)
    call muscle(int(four ),  I_limiter_type, first, last, tree)
    call muscle(int(five ),  I_limiter_type, first, last, tree)
    call muscle(int(six  ),  I_limiter_type, first, last, tree)
    
    ! ===> compute rhs
    call loop_on_quadtree_array(first, last, tree, compute_RHS)
    
    ! ==> update conservative variables
    call  loop_on_quadtree_array(first, last, tree, second_Ex_single)
    
    ! ===> update time
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