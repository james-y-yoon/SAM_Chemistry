MODULE cloudchem_Model

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Completely defines the model cloudchem
!    by using all the associated modules
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE cloudchem_Precision
  USE cloudchem_Parameters
  USE cloudchem_Global
  USE cloudchem_Function
  USE cloudchem_Integrator
  USE cloudchem_Rates
  USE cloudchem_Jacobian
  USE cloudchem_LinearAlgebra
  USE cloudchem_Monitor
  USE cloudchem_Util

END MODULE cloudchem_Model

