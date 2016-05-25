#!/bin/bash

# turn on logging and exit on error
set -e -x

sudo tee /etc/apt/sources.list.d/R.list <<EOF
# Change these lines if you don't want to use the main CRAN mirror.
# debian R upgrade
deb http://cran.r-project.org/bin/linux/debian squeeze-cran/
deb-src http://cran.r-project.org/bin/linux/debian squeeze-cran/

EOF


# ## THE FOLLOWING WAS NEEDED IN LENNY, BUT IT APPEARS NOT IN SQUEEZE
# sudo tee /etc/apt/preferences <<EOF
#  Package: *
#  Pin: release a=oldstable
#  Pin-Priority: 910

#  Package: *
#  Pin: release a=stable
#  Pin-Priority: 920

#  Package: *
#  Pin: release a=testing
#  Pin-Priority: 900

#  Package: *
#  Pin: release a=unstable
#  Pin-Priority: 800
# EOF

## test
# echo "force-confold" | sudo tee -a  /etc/dpkg/dpkg.cfg
# echo "force-confdef" | sudo tee -a  /etc/dpkg/dpkg.cfg

# add key to keyring so it doesn't complain
 gpg --keyserver pgp.mit.edu --recv-key 381BA480
 gpg -a --export 381BA480 > jranke_cran.asc
 sudo apt-key add jranke_cran.asc

# install the gfortran
sudo apt-get update
## when we went to Squeeze from lenny, 4.2 no longer installs. Is it needed? time will tell!! 
## sudo apt-get install --yes gfortran-4.2

# install R using the FRONTEND call to eliminate
# user interactive requests
sudo DEBIAN_FRONTEND=noninteractive apt-get install -t testing --yes --force-yes gcc
sudo DEBIAN_FRONTEND=noninteractive apt-get install -t testing --yes --force-yes r-base
sudo DEBIAN_FRONTEND=noninteractive apt-get install -t testing --yes --force-yes r-base-dev r-cran-hmisc

## rJava and latest Sun Java
## causing errors so dropping
# sudo DEBIAN_FRONTEND=noninteractive apt-get install -t unstable --yes --force-yes sun-java6-jdk sun-java6-jre r-cran-rjava

## get rJava working, by any means possible
# echo "### Hacked in to get rJava working ###" | sudo tee -a  /home/hadoop/.bashrc
# echo "export JAVA_HOME=/usr/lib/jvm/java-6-sun/jre" | sudo tee -a  /home/hadoop/.bashrc
# sudo env JAVA_HOME=/usr/lib/jvm/java-6-sun/jre R CMD javareconf

#install littler
sudo apt-get install -t testing littler

#some packages have trouble installing without this link
sudo ln -s /usr/lib/libgfortran.so.3 /usr/lib/libgfortran.so

# for the package update script to run
# the user hadoop needs to own the R library
sudo chown -R hadoop /usr/lib/R/library

