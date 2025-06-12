bitsadmin /transfer lanmsngrDownload https://github.com/lanmessenger/lanmessenger/releases/latest/download/lanmsngr-windows-x64.zip "%USERPROFILE%\Downloads\lanmsngr.zip"; `
Expand-Archive "$env:USERPROFILE\Downloads\lanmsngr.zip" "$env:USERPROFILE\LanMessenger" -Force; `
Start-Process "$env:USERPROFILE\LanMessenger\lanmsngr.exe"
