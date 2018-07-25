!=================================================================================================
!> CARTESIAN ADAPTIVE FIVE EQUATION MODEL
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/2D_CARFIVE
!================================================================================================= 
MODULE MODULE_MUSCL
    
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    use MODULE_GENERICMETHOD
    use MODULE_QUADTREE
    use MODULE_CFDUTILITY
    
!=========================== FACTORY PATTERN ===========================
    type, abstract  ::  limiter_funcs
    contains
    procedure(generic_funcs), deferred, pass(this)  :: l_funcs
    end type limiteR_funcs        
    
    abstract interface
        subroutine generic_funcs(this, a, b, omega, delta)
        use MODULE_PRECISION        
        import  :: limiter_funcs
        class(limiter_funcs), intent(in)    :: this
        real(rp), intent(in)                :: a, b, omega
        real(rp), intent(out)               :: delta
        end subroutine generic_funcs
    end interface
    
    type :: limiter_funcs_factory
    	private
        character(len = 20)             :: limiter_type
        class(limiter_funcs), pointer   :: limiter
    contains
    
    procedure   :: create_limiter
    
    end type limiter_funcs_factory   
    
    type, extends(limiter_funcs)    :: UST_super_bee
    contains
    procedure, pass(this)   :: l_funcs => UST_super_bee_limiter_funcs
    end type UST_super_bee     
!======================== END FACTORY PATTERN ==========================  
    
    contains
!==================================================================================================    
    function create_limiter(this, LIMITE) result(ptr)
    implicit none
    class(limiter_funcs_factory)    :: this
    class(limiter_funcs), pointer   :: ptr
    integer(ip), intent(in)  :: LIMITE
    
    if (LIMITE == 1) then 
        ! LIMITE = 1, Godunov's first order upwind method
        allocate(UST_super_bee :: this%limiter)
        this%limiter_type = "super bee"
        ptr => this%limiter
    else 
        write(*, *), "The program have not supported this limiter yet"
        write(*, *), "Stop."
        pause
        stop
    end if

    return
    end function create_limiter
!==================================================================================================
    subroutine UST_super_bee_limiter_funcs(this, a, b, omega, delta)
!
!   purpose: to compute a superbee type slope limiter delta
!
    implicit none
!
!   declaration of variables
!
    class(UST_super_bee),  intent(in)  :: this
    real(rp), intent(in)    :: a, b, omega
    real(rp), intent(out)   :: delta
    real(rp)  denor, phi, phir
!
      phi = (sign(one, a) + sign(one, b)) / two * max(min(omega*abs(a), abs(b)), min(abs(a), omega*abs(b)))
    delta = phi
!
    end subroutine UST_super_bee_limiter_funcs    
!==================================================================================================    
    subroutine muscle(iVal, limite, first, last, tree)
    implicit none

    integer(ip), intent(in)                                 :: iVal, limite
    integer(ip), intent(in)                                 :: first, last
    type(quadtree), dimension(first:last), intent(inout)    :: tree
    type (limiter_funcs_factory)                            :: factory
    class (limiter_funcs), pointer                          :: limiter => null()    
    real(rp)    :: dupw, dloc
    integer(ip) :: i, j, k
    
    real(rp), parameter :: tollim = tolerance
    real(rp), parameter :: omega = third
    
    real(rp)    :: delta, ratio
    real(rp)    :: ratio_inv, delta_inv
    
!> create limiter

    limiter => factory%create_limiter(limite)
    
    call loop_on_quadtree_array(first, last,  tree, MUSCL_single)
!    ! X_component
!    do j = 1, ny
!        do i = 2, nx-1
!            dupw = w(i  ,j,iVal) - w(i-1,j,iVal)
!            dloc = w(i+1,j,iVal) - w(i  ,j,iVal)
!            
!            if(abs(dupw).le.tollim) dupw = tollim*sign(one,dupw)
!            if(abs(dloc).le.tollim) dloc = tollim*sign(one,dloc)
!            
!            !delta = half*(one + omega)*dupw + half*(one - omega)*dloc      !(old version muscle center low accuracy, dont use any more)
!
!                ratio = dupw/dloc
!            ratio_inv = one/ratio
!            
!            ! Slope limiters used are:
!            ! 
!            ! LIMITE = 1, Godunov's first order upwind method           ( both are fird order, need tobe checked, not important, we don't use it any more)
!            ! LIMITE = 2, upwind second order method (non-monotone)     ( both are fird order, need tobe checked, not important, we don't use it any more)
!            ! LIMITE = 3, upwind TVD, with SUPERBEE type limiter        (old superbee riemann solver book)
!            ! LIMITE = 4, upwind TVD, with VAN LEER type limiter        (             riemann solver book)
!            ! LIMITE = 5, upwind TVD, with VAN ALBADA type limiter      (             riemann solver book)
!            ! LIMITE = 6, upwind TVD, with MINMOD type limiter  >>>>> NOT WORK ANY MORE <<<<<
!            ! LIMITE = 7, upwind TVD, with MINMAX type limiter  >>>>> NOT WORK ANY MORE <<<<<
!            ! LIMITE = 6, new superbee (high accuracy source: https://en.wikipedia.org/wiki/Flux_limiter)
!            call limiter%l_funcs(ratio, omega, delta )    
!            delta_inv = delta/ratio
!            
!            muscl_recons(i,j,ival)%x_l = w(i,j,ival) - fourth*((one + omega)*delta_inv*dupw + (one - omega)*delta*dloc)
!            muscl_recons(i,j,ival)%x_r = w(i,j,ival) + fourth*((one - omega)*delta_inv*dupw + (one + omega)*delta*dloc)
!
!        end do
!    end do

    ! y_component
!    do j = 2, ny-1
!        do i = 1, nx
!            dupw = w(i,j  ,iVal) - w(i,j-1,iVal)
!            dloc = w(i,j+1,iVal) - w(i,j  ,iVal)
!            
!            if(abs(dupw).le.tollim)dupw=tollim*sign(one,dupw)
!            if(abs(dloc).le.tollim)dloc=tollim*sign(one,dloc)
!            
!            !delta = half*(one + omega)*dupw + half*(one - omega)*dloc
!
!                ratio = dupw/dloc
!            ratio_inv = one/ratio
!            
!            ! Slope limiters used are:
!            ! 
!            ! LIMITE = 1, Godunov's first order upwind method
!            ! LIMITE = 2, upwind second order method (non-monotone)
!            ! LIMITE = 3, upwind TVD, with SUPERBEE type limiter
!            ! LIMITE = 4, upwind TVD, with VAN LEER type limiter
!            ! LIMITE = 5, upwind TVD, with VAN ALBADA type limiter
!            ! LIMITE = 6, upwind TVD, with MINMOD type limiter
!            ! LIMITE = 7, upwind TVD, with MINMAX type limiter
!            call limiter%l_funcs(ratio, omega, delta )    ; delta_inv = delta/ratio
!
!            muscl_recons(i,j,ival)%y_l = w(i,j,ival) - fourth*((one + omega)*delta_inv*dupw + (one - omega)*delta*dloc)
!            muscl_recons(i,j,ival)%y_r = w(i,j,ival) + fourth*((one - omega)*delta_inv*dupw + (one + omega)*delta*dloc)
            
!        end do
!    end do
    
!     muscl_recons(1,:,:) = muscl_recons(2,:,:)
!    muscl_recons(nx,:,:) = muscl_recons(nx-1,:,:)
!     muscl_recons(:,1,:) = muscl_recons(:,2,:)
!    muscl_recons(:,ny,:) = muscl_recons(:,ny-1,:)
    
!> delete limiter    
    deallocate(factory%limiter)
    return
    
! ===> CONTAINS <==================================================================================
    contains
! =================================================================================================
    
    subroutine  MUSCL_single(tree)
    implicit none
    type(quadtree), pointer, intent(inout)  :: tree
    real(rp)                                :: delta
    
    type(quadtree), pointer                 :: curr_C, adj_c
    real(rp), dimension(2)                  :: dr, w_curr, w_adj
    real(rp)                                :: omega, temp

    ! body
     omega = third              ! setup the parameter of muscl scheme
    curr_c => tree              ! current cell
    ! ===> north cell <============================================================================
    if ( associated(tree%adj_north))  then 
        adj_c => tree%adj_north    ! adjoint cell
    else
        adj_c => tree
    end if
    
    call muscl_single_edge(curr_c, adj_c)
    
    ! ===> east cell <=============================================================================
    if ( associated(tree%adj_east))  then 
        adj_c => tree%adj_east    ! adjoint cell
    else
        adj_c => tree
    end if
    
    call muscl_single_edge(curr_c, adj_c)
    
    ! ===> south <=================================================================================
    if ( associated(tree%adj_south))  then 
        adj_c => tree%adj_south    ! adjoint cell
    else
        adj_c => tree
    end if
    
    call muscl_single_edge(curr_c, adj_c)   
    
    ! ===> west <==================================================================================
    if ( associated(tree%adj_west))  then 
        adj_c => tree%adj_west    ! adjoint cell
    else
        adj_c => tree
    end if
    
    call muscl_single_edge(curr_c, adj_c)   
    
    return
    end subroutine MUSCL_single
    
    subroutine muscl_single_edge(curr_c, adj_c)
    implicit none
    type(quadtree), intent(inout)   :: curr_c, adj_c
    
    real(rp), dimension(2)          :: dr, w_curr
    real(rp)                        :: delta, temp
    
    ! body 
    delta = adj_c%data%w(iVal) - curr_c%data%w(iVal)
    
    dr(:) = half*(curr_c%pts(1)%coord(:) + curr_c%pts(2)%coord(:)) - curr_c%pts(5)%coord(:)
      
    ! ===> compute limiter funcs <=================================================================
    call weighted_least_square(iVal, curr_c, w_curr)
    call limiter%l_funcs(two*dot_product(w_curr,dr) - delta, delta, omega, temp)
    
    ! ===> compute muscl reconstruction <==========================================================
    curr_c%data%recons(iVal)%x_r = curr_c%data%w(iVal) + half*temp
    
    end subroutine muscl_single_edge
    
    end subroutine muscle
!==================================================================================================
!
      subroutine sbslic(r, omega, delta)
!
!     purpose: to compute a superbee type slope limiter delta
!
      implicit none
!
!     declaration of variables
!
      real(rp)  delta, denor, omega, phi, phir, r
!
      phi             = zero
      if(r.ge.zero)phi = two*r
      if(r.ge.half)phi = one
!     
      if(r.ge.one)then
         denor = one - omega + (one + omega)*r
         phir  = two/denor
         phi   = min(phir, r)
         phi   = min(phi, two)
      endif

!
      delta = phi!*delta
!
      end
!
!==================================================================================================
!
      subroutine superbee(r, omega, delta)
!
!     purpose: to compute a superbee type slope limiter delta
!
      implicit none
!z
!     declaration of variables
!
      real(rp)  delta, denor, omega, phi, phir, r
!
        phi = zero
        
        phi = max(zero, min(two*r, one), min(r, two))
      delta = phi
!
      end
!
!==================================================================================================     
      
!
      subroutine vlslic(r, omega, delta)
!
!     purpose: to compute a van leer type slope limiter delta
!
      implicit none
!
!     declaration of variables
!
      real(rp)  delta, denor, omega, phi, phir, r
!
      phi = zero
!
      if(r.ge.zero)then
         denor = one - omega + (one + omega)*r
         phir  = two/denor
         phi   = two*r/(one + r)
         phi   = min(phi, phir)
      endif
!
      delta    = phi!*delta
!
      end
!
!==================================================================================================
!
      subroutine vaslic(r, omega, delta)
!
!     purpose: to compute a van albada type slope limiter delta
!
      implicit none
!
!     declaration of variables
!
      real(rp)  delta, denor, omega, phi, phir, r
!
      phi = zero
!
      if(r.ge.zero)then
         denor = one - omega + (one + omega)*r
         phir  = two/denor
         phi   = r*(one + r)/(one + r*r)
         phi   = min(phi, phir)
      endif
!
      delta    = phi!*delta
!
      end
!
!==================================================================================================
!
      subroutine mislic(r, omega, delta)
!
!     purpose: to compute a minmod type slope limiter delta
!
      implicit none
!
!     declaration of variables
!
      real(rp)  delta, denor, omega, phi, phir, r
!
      phi             = zero
      if(r.ge.zero)phi = r
!
      if(r.ge.one)then
         denor = two*(one - omega + (one + omega)*r)
         phir  = four/denor
         phi   = min(one, phir)
      endif
!
      delta    = phi!*delta
!
      end
!
!==================================================================================================
!
      subroutine minmax(dupw, dloc, delta)
!
!     purpose: to compute a minmax type slope limiter delta.
!              this is the most diffusive of all limiters
!              for centred schemes
!
      implicit none
!
!     declaration of variables
!
      real(rp)  betal, betar, delta, dloc, dupw, signo
!
      betal = one
      betar = one
      signo = half*(sign(one,dupw) + sign(one,dloc))
      delta = signo*(min(betal*abs(dupw),betar*abs(dloc)))
!
      end
!
!================================================================================================== 
END MODULE MODULE_MUSCL
