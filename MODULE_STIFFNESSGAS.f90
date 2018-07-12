!=================================================================================================
!> CARTESIAN ADAPTIVE FIVE EQUATION MODEL
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/2D_CARFIVE
!================================================================================================= 
MODULE MODULE_STIFFNESSGAS
  
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    
    contains
!==================================================================================================
    function compute_gamma_mixture(a1, a2, gamma1, gamma2) result(res)
    implicit none
    real(rp), intent(in)    :: a1, a2, gamma1, gamma2
    real(rp)                :: res
    
    res = one/(a1/(gamma1 - one) + a2/(gamma2 - one)) + one
    return
    end function compute_gamma_mixture
!==================================================================================================  
    function compute_pi_mixture(a1, a2, pi1, pi2, gamma1, gamma2) result(res)
    implicit none
    real(rp), intent(in)    :: a1, a2, pi1, pi2, gamma1, gamma2
    real(rp)                :: res, gamma_mix
    
    gamma_mix = compute_gamma_mixture(a1, a2, gamma1, gamma2)
          res = (a1*pi1*gamma1/(gamma1 - one) + a2*pi2*gamma2/(gamma2 - one))*(gamma_mix - one)/gamma_mix
    return
    end function compute_pi_mixture
!================================================================================================== 
    function compute_rhoE(P, gamma_mix, pi_mix) result(res)
    implicit none
    real(rp), intent(in)    :: P, gamma_mix, pi_mix
    real(rp)                :: res
    
    res = (P + gamma_mix*pi_mix)/(gamma_mix - one)
    return
    end function compute_rhoE
!================================================================================================== 
    function compute_P(gamma_mix, rho_e, pi_mix) result(res)
    implicit none
    real(rp), intent(in)    :: gamma_mix, rho_e, pi_mix
    real(rp)                :: res
    
    res = (gamma_mix - one)*rho_e - gamma_mix*pi_mix
    return
    end function compute_P
!================================================================================================== 
    function compute_c_mixture(gamma_mix, rho, p, pi_mix) result(res)
    implicit none
    real(rp), intent(in)    :: gamma_mix, rho, p, pi_mix
    real(rp)                :: res
    
    res = sqrt(gamma_mix*(p + pi_mix)/rho)
    return
    end function compute_c_mixture
!==================================================================================================    
END MODULE MODULE_STIFFNESSGAS    