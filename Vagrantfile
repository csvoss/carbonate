$script = <<SCRIPT
sudo apt-get update
sudo apt-get install -y git curl python-pip 
cd /tmp
pip install Django South django-debug-toolbar django-extensions
SCRIPT

VAGRANTFILE_API_VERSION = "2"
Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  config.vm.box = "hashicorp/precise32"
  config.vm.provision :shell, :inline => $script
  config.vm.network :forwarded_port, guest: 4000, host: 4000
  config.vm.provider "virtualbox" do |v|
      v.customize ["setextradata", :id, "VBoxInternal2/SharedFoldersEnableSymlinksCreate/v-root", "1"]
  end
end
