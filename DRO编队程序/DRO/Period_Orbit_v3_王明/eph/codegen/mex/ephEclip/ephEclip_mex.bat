@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2016a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2016a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=ephEclip_mex
set MEX_NAME=ephEclip_mex
set MEX_EXT=.mexw64
call "C:\PROGRA~1\MATLAB\R2016a\sys\lcc64\lcc64\mex\lcc64opts.bat"
echo # Make settings for ephEclip > ephEclip_mex.mki
echo COMPILER=%COMPILER%>> ephEclip_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> ephEclip_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> ephEclip_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> ephEclip_mex.mki
echo LINKER=%LINKER%>> ephEclip_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> ephEclip_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> ephEclip_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> ephEclip_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> ephEclip_mex.mki
echo BORLAND=%BORLAND%>> ephEclip_mex.mki
echo OMPFLAGS= >> ephEclip_mex.mki
echo OMPLINKFLAGS= >> ephEclip_mex.mki
echo EMC_COMPILER=lcc64>> ephEclip_mex.mki
echo EMC_CONFIG=optim>> ephEclip_mex.mki
"C:\Program Files\MATLAB\R2016a\bin\win64\gmake" -B -f ephEclip_mex.mk
