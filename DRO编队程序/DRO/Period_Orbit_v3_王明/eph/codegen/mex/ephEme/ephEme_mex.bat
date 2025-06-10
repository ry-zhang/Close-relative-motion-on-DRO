@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2016a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2016a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=ephEme_mex
set MEX_NAME=ephEme_mex
set MEX_EXT=.mexw64
call "C:\PROGRA~1\MATLAB\R2016a\sys\lcc64\lcc64\mex\lcc64opts.bat"
echo # Make settings for ephEme > ephEme_mex.mki
echo COMPILER=%COMPILER%>> ephEme_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> ephEme_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> ephEme_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> ephEme_mex.mki
echo LINKER=%LINKER%>> ephEme_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> ephEme_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> ephEme_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> ephEme_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> ephEme_mex.mki
echo BORLAND=%BORLAND%>> ephEme_mex.mki
echo OMPFLAGS= >> ephEme_mex.mki
echo OMPLINKFLAGS= >> ephEme_mex.mki
echo EMC_COMPILER=lcc64>> ephEme_mex.mki
echo EMC_CONFIG=optim>> ephEme_mex.mki
"C:\Program Files\MATLAB\R2016a\bin\win64\gmake" -B -f ephEme_mex.mk
