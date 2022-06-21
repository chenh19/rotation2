#!/bin/bash

# Canon Color ImageCLASS MF735Cdw driver for linux
cd ~/Downloads
wget -q -O linux-UFRII-drv-v550-us-00.tar.gz https://pdisp01.c-wss.com/gdl/WWUFORedirectTarget.do?id=MDEwMDAwOTIzNjEx&cmp=ABR&lang=EN && sleep 10
tar -xf linux-UFRII-drv-v550-us-00.tar.gz && sleep 1 && rm linux-UFRII-drv-v550-us-00.tar.gz
sudo bash ./linux-UFRII-drv-v550-us/install.sh && sleep 1 && rm -rf ./linux-UFRII-drv-v550-us/
sudo apt-get update && sudo apt install -f -y && sudo apt-get autoremove -y && sudo apt-get clean
