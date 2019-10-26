# FILESTORE NFS SETUP NOTES

## Create Filestore instance

`gcloud filestore instances create finngen-nfs --project phewas-development --zone europe-west1-b --tier=PREMIUM --file-share=name="vol1",capacity=4TB --network=name="default"`

To mount on a VM in the same project:

```
sudo mkdir /mnt/nfs
sudo mount 10.179.247.250:/vol1 /mnt/nfs/
```

where the IP address is the address of the Filestore instance.

## Set up a VPN server to be able to mount the file system locally

Following/applying [this](https://medium.com/teendevs/setting-up-an-openvpn-server-on-google-compute-engine-9ff760d775d9)

We're using port 1194 and UDP instead of 443 TCP in the instructions. Make sure you have a firewall rule allowing UDP connections to port 1194 on the VPN VM.

Set up VPN server on the VM:

```
sudo apt-get update
sudo apt-get install openvpn easy-rsa emacs -y
sudo make-cadir /opt/openvpn-ca
sudo chown -R jkarjala /opt/openvpn-ca/
cd /opt/openvpn-ca/
emacs vars
```

Set e.g.

```
export KEY_COUNTRY="FI"
export KEY_PROVINCE="HEL"
export KEY_CITY="Helsinki"
export KEY_ORG="FinnGen"
export KEY_EMAIL="humgen-servicedesk@helsinki.fi"
export KEY_OU="FIMM"
export KEY_NAME="server"
```

Then

```
ln -s openssl-1.0.0.cnf openssl.cnf
source vars
./clean-all
./build-ca
./build-key-server server
./build-dh
openvpn --genkey --secret keys/tiv.key

./build-key client_jk
./build-key client_mk

cd keys
sudo cp ca.crt server.crt server.key tiv.key dh2048.pem /etc/openvpn/
gunzip -c /usr/share/doc/openvpn/examples/sample-config-files/server.conf.gz | sudo tee /etc/openvpn/server.conf
sudo emacs /etc/openvpn/server.conf
```

Make sure the following are set:

```
port 1194

;proto tcp
proto udp

tls-auth tiv.key 0 # This file is secret                                                                                                                                                                                                      
key-direction 0

cipher AES-256-CBC
auth SHA256

user nobody
group nogroup

push "redirect-gateway def1 bypass-dhcp"

push "dhcp-option DNS 208.67.222.222"
push "dhcp-option DNS 208.67.220.220"
```

```
sudo emacs /etc/sysctl.conf
```

Set

```
net.ipv4.ip_forward=1
```

Then

```
sudo sysctl -p
ip route | grep default
#default via 10.132.0.1 dev ens4 proto dhcp src 10.132.0.63 metric 100

sudo emacs /etc/ufw/before.rules
```

Add

```
# OPENVPN
# NAT Table
*nat
:POSTROUTING ACCEPT [0:0] 
# OpenVPN client traffic
-A POSTROUTING -s 10.8.0.0/8 -o ens4 -j MASQUERADE
COMMIT
# OPENVPN
```

Not mentioned in the instructions, this is also required:

```
sudo iptables -t nat -A POSTROUTING -s 10.8.0.0/8 -o ens4 -j MASQUERADE
```

```
sudo emacs /etc/default/ufw
```

Set

```
DEFAULT_FORWARD_POLICY="ACCEPT"
```

Finally start the server and have it start automatically on startup:

```
sudo systemctl start openvpn@server
sudo systemctl status openvpn@server
sudo systemctl enable openvpn@server
```

## Create client configs

```
sudo mkdir -p /opt/clients/files
sudo chown -R jkarjala /opt/clients
chmod 700 /opt/clients/files
cp /usr/share/doc/openvpn/examples/sample-config-files/client.conf /opt/clients/base.conf
emacs /opt/clients/base.conf
```

Make sure the following are set (34.76.184.95 is the IP of the VM):

```
;proto tcp
proto udp

remote 34.76.184.95 1194

user nobody
group nogroup

# ca ca.crt                                                                                                                                                                                                                                   
# cert client.crt                                                                                                                                                                                                                             
# key client.key                                                                                                                                                                                                                              

cipher AES-256-CBC
auth SHA256

key-direction 1
```

```
emacs /opt/clients/gen_config.sh
```

Paste this

```
#!/bin/bash                                                                                                                                                                                                                                   
KEY_DIR=/opt/openvpn-ca/keys
OUTPUT_DIR=/opt/clients/files
BASE_CONFIG=/opt/clients/base.conf
cat ${BASE_CONFIG} \
    <(echo -e '<ca>') \
    ${KEY_DIR}/ca.crt \
    <(echo -e '</ca>\n<cert>') \
    ${KEY_DIR}/${1}.crt \
    <(echo -e '</cert>\n<key>') \
    ${KEY_DIR}/${1}.key \
    <(echo -e '</key>\n<tls-auth>') \
    ${KEY_DIR}/tiv.key \
    <(echo -e '</tls-auth>') \
    > ${OUTPUT_DIR}/${1}.ovpn
```

Then (change your_initials to your initials):

```
chmod 700 /opt/clients/gen_config.sh

cd /opt/clients/
./gen_config.sh client_your_initials

emacs files/client_your_initials.ovpn
```

Set

```
tls-auth tiv.key 1
```

## Get the configs to your local machine

Change your_initials to your initials

```
mkdir -p ~/.vpn
cd ~/.vpn
gcloud compute scp vpn:/opt/clients/files/client_your_initials.ovpn . --project phewas-development
gcloud compute scp vpn:/opt/openvpn-ca/keys/tiv.key . --project phewas-development
```

Then drag the .ovpn file to Tunnelblick and good to go.

On MacOS:

```
sudo mkdir -p /mnt/nfs
sudo mount_nfs -o resvport 10.179.247.250:/vol1 /mnt/nfs/
```

If you at some point can't connect to the VPN or can't SSH to the VPN server VM, your routing table may need to be flushed:

```
sudo ifconfig en0 down
sudo route flush
sudo ifconfig en0 up
```
