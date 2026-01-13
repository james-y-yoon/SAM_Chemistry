!WRF:MODEL_LAYER:CONSTANTS
!
!bloss: Trim WRFv4.2.2 version of this module to contain only the needed constants.

 MODULE module_model_constants

   REAL    , PARAMETER :: RE_QC_BG     = 2.49E-6     ! effective radius of cloud for background (m)
   REAL    , PARAMETER :: RE_QI_BG     = 4.99E-6     ! effective radius of ice for background (m)
   REAL    , PARAMETER :: RE_QS_BG     = 9.99E-6     ! effective radius of snow for background (m)

 END MODULE module_model_constants
