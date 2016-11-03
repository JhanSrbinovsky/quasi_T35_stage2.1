module cable_decs_mod

  type CABLE_file
    logical :: want = .false. 
    integer :: funit
    character(len=199) :: folder
    character(len=190) :: filename
    logical :: vars = .false.
    integer :: itau=0,ftau=0,mtau=0 
  End type CABLE_file
  
End module cable_decs_mod
