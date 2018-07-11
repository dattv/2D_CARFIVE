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
    
END MODULE MODULE_CFDMAINDATA    