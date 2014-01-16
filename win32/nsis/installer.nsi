
!include "MUI2.nsh"
!addincludedir "../win32/nsis/include"

Name "pyseq"
OutFile "pyseq-installer-test.exe"
InstallDir "$PROGRAMFILES\pyseq"
InstallDirRegKey HKCU "Software\pyseq" ""

RequestExecutionLevel user

!define MUI_ABORTWARNING

;!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES

!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES

Section "Pyseq" SecDummy
  SetOutPath "$INSTDIR"
  !include "install.files"
  File simple.exe
  File /r tcl
  WriteRegStr HKCU "Software\pyseq" "" $INSTDIR
  WriteUninstaller "$INSTDIR\Uninstall.exe"
SectionEnd

;Language strings
LangString DESC_SecDummy ${LANG_ENGLISH} "A test section."

;Assign language strings to sections
!insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
!insertmacro MUI_DESCRIPTION_TEXT ${SecDummy} $(DESC_SecDummy)
!insertmacro MUI_FUNCTION_DESCRIPTION_END

Section "Uninstall"
  !include "uninstall.files"
  Delete "$INSTDIR\simple.exe"
  RMDir /r "$INSTDIR\tcl"
  Delete "$INSTDIR\Uninstall.exe"
  RMDir "$INSTDIR"
  DeleteRegKey /ifempty HKCU "Software\pyseq"
SectionEnd

Function .onInit
  ;TODO check whether already installed
FunctionEnd

