  !-----------------------------------------------------------------------
  !
  ! CUSTOM MODELS FOR SPECFEM3D_GLOBE
  ! Add to the case selection in: src/shared/get_model_parameters.F90
  !
  !
  ! custom_1d_transversely_isotropic_prem: 1D TISO PREM that does not honor
  ! PREM moho and allows for 3D crustal variations, allows 1D TISO PREM to be
  ! used for GLL inversions, otherwise mesh variations will cause nelement
  ! mismatches when calling model='gll'
  !
  ! custom_gll: to be used with `custom_1d_transversely_isotropic_prem`, turns
  ! off some gll hardcoding and makes the mesh look the same as the custom 1D
  ! TISO prem so that using model='custom_gll' allows model updates in SeisFlows
  !
  ! Important Parameters:
  ! CASE_3D: If True, stretches elements in the crustal layers to get more GLL 
  !     points per km in the upper crust to allow for 3d heterogeneities in the 
  !     upper crust. 
  !
  ! CRUSTAL: If True, will superimpose a 3D crustal velocity model onto the 
  !     crustal elements. This happens by extending 1D mantle values into the 
  !     crust adn then overprinting with a 3D crustal model. If False then a 1D 
  !     crustal model will be taken from the 1D reference model
  !
  ! ONE_CRUST: To increase stability and allow for cheaper simulations, 1D 
  !     models can average the crust into a single layer (True) instead of 2 
  !     (False). 
  !     WARNING: Be extremely cautious when mixing this with custom
  !     Moho depths etc., I was getting weird artefacting when ONE_CRUST=True
  !     where parts of the crust were randomly one or two layer. 
  !     
  !
  !-----------------------------------------------------------------------
  case ('custom_1d_transversely_isotropic_prem')
    CASE_3D = .false.    ! true = more GLL points / km in upper crust
    CRUSTAL = .false.    ! true: impose 3D crustal model; false: 1D crust
    ONE_CRUST = .true.   ! true = 1 crust layer; false = 2 crust layers

    TRANSVERSE_ISOTROPY = .true.  ! true = vpv, vph, vsv, vsh, rho
    HONOR_1D_SPHERICAL_MOHO = .true.

    ! These are default values assigned prior to model selection and have not
    ! been changed, but are left here as reference and incase something changes
    ! in future SPECFEM verisons
    REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
    HETEROGEN_3D_MANTLE = .false.
    MODEL_3D_MANTLE_PERTUBATIONS = .false.
    MODEL_GLL = .false.
    MODEL_GLL_TYPE = 0
    MODEL_3D_MANTLE_PERTUBATIONS = .false.  
    ATTENUATION_3D = .false.
    CEM_REQUEST = .false.
    CEM_ACCEPT  = .false.
    EMC_MODEL = .false.

  case ('custom_gll')
    ! Make sure these match `custom_1d_transversely_isotropic_prem`
    CASE_3D = .false.    ! default .true.
    CRUSTAL = .false.    ! default .true.
    ONE_CRUST = .true.   ! default .true.

    MODEL_3D_MANTLE_PERTUBATIONS = .false.   ! default .true.

    ! The following parameters have not been modified
    TRANSVERSE_ISOTROPY = .true.  ! same as reference model

    MODEL_GLL = .true.
    MODEL_GLL_TYPE = 2 ! (2 == tiso) input model files are tiso

    REFERENCE_1D_MODEL = GLL_REFERENCE_1D_MODEL
    HETEROGEN_3D_MANTLE = .false.
    MODEL_3D_MANTLE_PERTUBATIONS = .false.
    MODEL_GLL = .false.
    MODEL_GLL_TYPE = 0
    ATTENUATION_3D = .false.
    CEM_REQUEST = .false.
    CEM_ACCEPT  = .false.
    EMC_MODEL = .false.
