@echo off
call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\VsDevCmd.bat"

REM q:quiet m:minimum n:normal
set "PARAMS=/m /nr:false /nologo /verbosity:q /t:Build /p:Platform=x86"

rem ========== build cmds ==========
call "cleanall.bat"
call :build siglibs.sln Release
call "clean.bat"
pause
goto end

rem ========== functions ==========
:build
SETLOCAL
@echo Building: %1 %2...
MSBuild %1 %PARAMS% /p:Configuration=%2
ENDLOCAL
goto end

:end
