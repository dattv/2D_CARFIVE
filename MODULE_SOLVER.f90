!=================================================================================================
!> CARTESIAN ADAPTIVE FIVE EQUATION MODEL
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/2D_CARFIVE
!================================================================================================= 
MODULE MODULE_SOLVER
    
    use MODULE_PRECISION
    
    type, abstract  :: abstract_solver
    contains
    procedure(generic_solver), deferred, pass(this) :: p_solver
    end type abstract_solver
    
    abstract interface
        subroutine generic_solver(this)
        use MODULE_PRECISION
        import  :: abstract_solver
        class(abstract_solver), intent(in)   :: this
        end subroutine generic_solver
    end interface

    type :: solver_factory
    	private
        character(len = 20)             :: c_solver_type
        integer(ip)                     :: i_solver_type
        class(abstract_solver), pointer :: 
    end type solver_factory
END MODULE MODULE_SOLVER    