module chemistry_string_limits
   implicit none
   private

!  CHEMKIN/NASA-7 thermo cards reserve columns 1:18 for the species
!  identifier. Mechanism and auxiliary lines are modern free-form text,
!  so keep their parsing workspace large enough for long reaction names.
   integer, parameter, public :: species_name_len  = 18
   integer, parameter, public :: mechanism_line_len = 256

end module chemistry_string_limits
