# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.

Vagrant.configure(2) do |config|

# Every Vagrant development environment requires a box. You can search for
# boxes at https://atlas.hashicorp.com/search.

config.vm.box = "scientific-linux-6"

# The url from where the 'config.vm.box' box will be fetched if it
# doesn't already exist on the user's system.

config.vm.box_url = "http://lyte.id.au/vagrant/sl6-64-lyte.box"

# Create a forwarded port mapping which allows access to a specific port
# within the machine from a port on the host machine. In the example below,
# accessing "localhost:8080" will access port 80 on the guest machine.
# config.vm.network "forwarded_port", guest: 80, host: 8080

# Create a private network, which allows host-only access to the machine
# using a specific IP.
# config.vm.network "private_network", ip: "192.168.33.10"

# Create a public network, which generally matched to bridged network.
# Bridged networks make the machine appear as another physical device on
# your network.
# config.vm.network "public_network"

# Share an additional folder to the guest VM. The first argument is
# the path on the host to the actual folder. The second argument is
# the path on the guest to mount the folder. And the optional third
# argument is a set of non-required options.

# config.vm.synced_folder "/Users/colemanderr/repo/test-sci-linux-vm/required/", "./required/"

# Provider-specific configuration so you can fine-tune various
# backing providers for Vagrant. These expose provider-specific options.
# Example for VirtualBox:
#
# config.vm.provider "virtualbox" do |vb|
#   # Display the VirtualBox GUI when booting the machine
#   vb.gui = true
#
#   # Customize the amount of memory on the VM:
#   vb.memory = "1024"
# end
#
# View the documentation for the provider you are using for more
# information on available options.

# Define a Vagrant Push strategy for pushing to Atlas. Other push strategies
# such as FTP and Heroku are also available. See the documentation at
# https://docs.vagrantup.com/v2/push/atlas.html for more information.
# config.push.define "atlas" do |push|
#   push.app = "YOUR_ATLAS_USERNAME/YOUR_APPLICATION_NAME"
# end

# Enable provisioning with a shell script. Additional provisioners such as
# Puppet, Chef, Ansible, Salt, and Docker are also available. Please see the
# documentation for more information about their specific syntax and use.
 


#config.vm.provision "file", source: "./vagrant_bash_profile", destination: ".bash_profile"

config.vm.provision "file", source: "./required/get-pip.py", destination: "get-pip.py"

config.vm.provision "shell", inline: <<-SHELL

source .bash_profile
mkdir bin

# adds nano as editor
# echo "Installing nano editor"
sudo yum install -y nano.x86_64

# adds devel version of python necessary for gcc compile of cutadapt
# sudo yum install -y  python-devel.x86_64

# curl -L http://install.perlbrew.pl | bash
# to upgrade to perl 5.18.4
# source ~/perl5/perlbrew/etc/bashrc
# to make perlbrew executable
# perlbrew install perl-5.18.4
# Adds latest version of perl, takes a while.
# perlbrew switch perl-5.18.4

# python get-pip.py    
# pip install cutadapt

# installing fastx_toolkit : instructions from http://hannonlab.cshl.edu/fastx_toolkit/install_centos.txt
# sudo yum install pkgconfig.x86_64 gcc.x86_64 gcc-c++.x86_64 wget.x86_64
# wget http://cancan.cshl.edu/labmembers/gordon/files/libgtextutils-0.6.tar.bz2
# tar -xjf libgtextutils-0.6.tar.bz2
# cd libgtextutils-0.6
# ./configure
# make
# sudo make install
# cd ..
# rm libgtextutils-0.6.tar.bz2

# export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH
# wget http://cancan.cshl.edu/labmembers/gordon/files/fastx_toolkit-0.0.12.tar.bz2 
# tar -xjf fastx_toolkit-0.0.12.tar.bz2 
# cd fastx_toolkit-0.0.12
# ./configure
# make
# sudo make install
# cd ..
# rm fastx_toolkit-0.0.12.tar.bz2

#getting devel version of zlib for installing FLASH
sudo yum install -y zlib-devel.x86_64

#getting flash 1.2.11
wget http://sourceforge.net/projects/flashpage/files/FLASH-1.2.11.tar.gz
tar xvfz FLASH-1.2.11.tar.gz
cd FLASH-1.2.11
make
cd ..
cp FLASH-1.2.11/flash bin/
rm FLASH-1.2.11.tar.gz


SHELL

end
