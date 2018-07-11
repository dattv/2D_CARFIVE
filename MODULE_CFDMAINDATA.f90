!=================================================================================================
!> CARTESIAN ADAPTIVE FIVE EQUATION MODEL
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/2D_CARFIVE
!=================================================================================================  
MODULE MODULE_CFDMAINDATA
    
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    
    ! ===> CLASS PROPERTIES <======================================================================
    type :: properties
        real(rp)    ::                    gamma = 1.4_rp        ;
	    real(rp)    ::                        b = 1._rp         ;
	    real(rp)    ::                      rho = 1.3_rp        ;
	    real(rp)    ::           r_gas_constant = 8.3144598_rp  ;
	    real(rp)    ::                       mu = 18.27_rp      ;
	    real(rp)    ::                   lambda = 1.512041288_rp;
	    real(rp)    ::   c_sutherland_constant  = 120._rp       ;
	    real(rp)    ::              temperature = 291.15_rp     ;
        real(rp)    ::                       pi = MPI           ;
        real(rp)    ::                       cp = 1._rp         ;
    end type properties
    ! ===> END CLASS PROPERTIES <==================================================================
    
    integer(ip) :: time_step_max    ! MAXIMUM PHYSICAL TIME STEPS
    real(rp)    :: CFL              ! CORRANT FEDRIC LEVI NUMBER
    real(rp)    :: t_final          ! FINAL TIME
    real(rp)    :: time             !
    real(rp)    :: dt               !

    ! ===> REFERENCE QUANTITIES <==================================================================
    real(rp)    :: M_inf    , &
                   Rho_inf  , &
                   U_inf    , &
                   P_inf    , &
                   L_inf
    
    ! ===> RATIO FO SPECIFIC HEAT (1.4 FRP AIR) <==================================================
    real(rp)    :: gamma = 1.4_rp
    
    character(len = 80) :: C_inviscid_flux
    integer(ip)         :: I_inviscid_flux
    
    character(len = 80) :: C_limiter_type
    integer(ip)         :: I_limiter_type
    
    integer(ip)         :: frequency_dump
    
    type(properties), dimension(:), pointer :: matInfo              
    
END MODULE MODULE_CFDMAINDATA    