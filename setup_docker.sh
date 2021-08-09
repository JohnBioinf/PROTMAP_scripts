#!/bin/bash

set -e

apt-get update
apt-get install -y unzip
apt-get install -y inkscape

# Specific for docker! dbus should run on any other distro out of box
dbus-uuidgen > /var/lib/dbus/machine-id
mkdir -p /var/run/dbus
dbus-daemon --config-file=/usr/share/dbus-1/system.conf --print-address

bin="$HOME/bin"
if [ ! -d "$bin" ]; then
	mkdir "$bin"
fi

# conda
conda config --append channels conda-forge
conda config --append channels bioconda
conda config --append channels r
conda env create -f conda_protmap.yml

# comet
wget -nv -P "$bin" https://sourceforge.net/projects/comet-ms/files/comet_2019014.zip
unzip -d "$bin" "$bin/comet_2019014.zip"
chmod +x "$bin/comet.2019014.linux.exe"
rm "$bin/comet_2019014.zip"

# chrome
# install normal chrome, to make sure all depencies are there

apt-get update 
apt-get install -y gnupg 
wget -q -O - https://dl-ssl.google.com/linux/linux_signing_key.pub | apt-key add - 
sh -c 'echo "deb [arch=amd64] http://dl.google.com/linux/chrome/deb/ stable main" >> /etc/apt/sources.list.d/google.list' 
apt-get update 
apt-get install -y google-chrome-stable fonts-ipafont-gothic fonts-wqy-zenhei fonts-thai-tlwg fonts-kacst fonts-freefont-ttf libxss1 --no-install-recommends
rm -rf /var/lib/apt/lists/*

# specific chrome and chromedriver version to make sure everthing runs correct
wget -nv -O "$bin/chrome-linux.zip" 'https://www.googleapis.com/download/storage/v1/b/chromium-browser-snapshots/o/Linux_x64%2F827102%2Fchrome-linux.zip?generation=1605233458736188&alt=media'
unzip -d "$bin" "$bin/chrome-linux.zip"
rm "$bin/chrome-linux.zip"

wget -nv -O "$bin/chromedriver_linux64.zip" 'https://www.googleapis.com/download/storage/v1/b/chromium-browser-snapshots/o/Linux_x64%2F827102%2Fchromedriver_linux64.zip?generation=1605233463367665&alt=media'
unzip -d "$bin" "$bin/chromedriver_linux64.zip"
cp "$bin/chromedriver_linux64/chromedriver" "$bin"
rm "$bin/chromedriver_linux64.zip"
