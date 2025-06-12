Set-ExecutionPolicy Bypass -Scope Process -Force; `
[System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072; `
iex ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))






choco install lanmessenger -y

choco uninstall lanmessenger -y


posle WinKey + R -> pisi %appdata% izbrisi sve od Lan ms app
