  Program FBC

  Use ctrl

  Implicit None 

  Call Initialize()
  if (.not.(flag_umbrella)) then
    Call simulation()
  else
    Call umbrella()
  end if
  
  End Program
