#include "MAPL.h"

MODULE fv_shared

  USE ESMF
  USE MAPL, ONLY : MAPL_FieldGet, MAPL_Verify

  IMPLICIT NONE
  !
  !  Functions used in both AdvCore_GridCompMod and DynCore_GridCompMod
  !
CONTAINS

     function get_short_name(field, rc) result(short_name)
      type(ESMF_Field) :: field
      integer, intent(out) :: rc
      character(len=:), allocatable :: short_name

      character(len=:), allocatable :: standard_name
      integer :: status

      call MAPL_FieldGet(field, standard_name=standard_name, _RC)
      select case (trim(standard_name))
      case ("specific_humidity")
         short_name = "Q"
      case ("mass_fraction_of_large_scale_cloud_liquid_water")
         short_name = "QLLS"
      case ("mass_fraction_of_convective_cloud_liquid_water")
         short_name = "QLCN"
      case ("mass_fraction_of_large_scale_cloud_ice_water")
         short_name = "QILS"
      case ("mass_fraction_of_convective_cloud_ice_water")
         short_name = "QICN"
      case ("large_scale_cloud_area_fraction")
         short_name = "CLLS"
      case ("convective_cloud_area_fraction")
         short_name = "CLCN"
      case ("mass_fraction_of_rain")
         short_name = "QRAIN"
      case ("mass_fraction_of_snow")
         short_name = "QSNOW"
      case ("mass_fraction_of_graupel")
         short_name = "QGRAUPEL"
      case default
         ! _FAIL("Unrecognized standard_name: " // trim(standard_name))
         short_name = trim(standard_name)
      end select

      _RETURN(_SUCCESS)
   end function get_short_name

   function field_is_cloud_water_species(field_name) result(is_cloud_water_species)
      character(len=*), intent(in) :: field_name
      logical :: is_cloud_water_species

      is_cloud_water_species = .false.
      if ( &
           (trim(field_name) == "Q") .or. &
           (trim(field_name) == "QLCN") .or. &
           (trim(field_name) == "QLLS") .or. &
           (trim(field_name) == "QICN") .or. &
           (trim(field_name) == "QILS") .or. &
           (trim(field_name) == "CLCN") .or. &
           (trim(field_name) == "CLLS") .or. &
           (trim(field_name) == "NCPL") .or. &
           (trim(field_name) == "NCPI") .or. &
           (trim(field_name) == "QRAIN") .or. &
           (trim(field_name) == "QSNOW") .or. &
           (trim(field_name) == "QGRAUPEL")) then
         is_cloud_water_species = .true.
      end if
   end function field_is_cloud_water_species

   function is_name_in_list(name, list) result(is_in_list)
      character(len=*), intent(in) :: name
      character(len=ESMF_MAXSTR), intent(in) :: list(:)
      logical :: is_in_list

      integer :: n

      is_in_list = .false.
      do n = 1, size(list)
         if (trim(name) == trim(list(n))) then
            is_in_list = .true.
            exit
         end if
      end do
   end function is_name_in_list

END MODULE fv_shared
