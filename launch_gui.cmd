@echo off
setlocal
pushd "%~dp0"
if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" -X utf8 -m hischeck.gui %*
) else (
    python -X utf8 -m hischeck.gui %*
)
if errorlevel 1 pause
popd
