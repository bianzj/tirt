; Extra installer steps for TiRT: a Start-menu shortcut to the PDF documents.
!macro customInstall
  CreateShortCut "$SMPROGRAMS\TiRT 文档.lnk" "$INSTDIR\docs"
!macroend

!macro customUnInstall
  Delete "$SMPROGRAMS\TiRT 文档.lnk"
!macroend