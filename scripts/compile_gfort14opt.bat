@rem     ######################################################
title    ### Compiling by gfortran (tdm-gcc-64): opt-mode   ###
@rem     ######################################################
rem      gfortran.exe�E�̏ꏊ
cd /d %~dp0
set      GNUHOME=C:\gfortran-14.2.0-msvcrt\mingw64
set      MPIHOME="C:\MSMPI\MSMPI64\Microsoft HPC Pack 2008 SDK"

set      SCLIB=..\libSpeedCHEM64.a

if       exist %SCLIB% (del %SCLIB% *.mod)

set      SRCS=working_precision.f90 ^
           SCutilities.f ^
           SCsparse_definitions.f90 ^
           SCsparse.f ^
           sparse_MPI.f ^
           dvode_f90_m.f90 ^
           SCmodule.f ^
           chemkin_module.f90 ^
           SCconV.f90 ^
           SCsetup.f ^
           SCallocate.f ^
           SCcklink.f90 ^
           gam.f ^
           gamsub.f ^
           opkdmain.f ^
           opkda1.f ^
           opkda2.f ^
           ddaspk.f ^
           rodas.f ^
           vode.f ^
           MEBDFSO.f ^
           rowmap.f ^
           radau5.f ^
           radaua.f ^
           chemistry_input.f90 ^
           SCbroadcast.f90

@rem     ######################################################
set      PATH=%GNUHOME%\bin;%PATH%
set      MPI_INC=%MPIHOME%\Include

set      FCFLAG=-O3 -fconvert=big-endian -I%MPI_INC% -m64 ^
               -fdefault-double-8 -fdefault-real-8 ^
               -ffixed-line-length-none -malign-double -fallow-argument-mismatch
@rem     ######################################################
@echo    Compiling now to get %SCLIB% ...
gfortran %FCFLAG% -c %SRCS% > compile_gft.log 2>&1
ar cr    %SCLIB% *.o >> compile_gft.log 2>&1
ranlib   %SCLIB% >> compile_gft.log 2>&1
@rem     ######################################################
@if      exist %SCLIB% (echo    Compiling was successfully done !!! )
del      *.o
pause