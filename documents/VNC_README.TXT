

================================================================
# TO CREATE OPENSTACK SERVER WITH VNC SERVER INSTALLED:

# Here's how to create an OpenStack virtual machine with VNC Server already installed  
# and ready to use, so that you can log into it remotely and use its graphical user 
# interface, just like a normal computer:

# First, make sure the custom CentOS7_GUI_VNCServer image (which I created) is 
# available in OpenStack:
openstack image list

# You might also want to look at the flavor list to help you choose a flavor (for
# how powerful your virtual machine will be):
openstack flavor list

# Then, create a new openstack server from that image, just like you normally would, 
# but you have to make sure you also include the security group: vnc_security_group. 
# Choose your openstack flavour. Replace "username" (after --key-name) with your own 
# username. Replace "your_server_name" with whatever name you choose:

openstack server create \
	--image CentOS7_GUI_VNCServer \
  --security-group default \
  --security-group vnc_security_group \
	--key-name username \
	--nic net-id=Private \
	--flavor mm1.medium \
	your_server_name

# Wait a bit and then check to make sure your server was created successfully:
openstack server list

# Next you need to assign a floating ip to your server.
# Look to see if there are any free floating ips:
openstack floating ip list
# If there are none available (i.e. that don't haev a fixed ip address assigned to it), get a new one: 
openstack floating ip create External
# assign a free floating ip to your server (replace your_server_name and 10.20.0.XXX with your specifics)
openstack server add floating ip your_server_name 10.20.0.XXX

================================================================
# TO CONNECT TO YOUR OPENSTACK VNC SERVER FROM YOUR OWN OFFICE COMPUTER:


----------
1. Download the client VNC Viewer.

The VNC client we're going to use is TigerVNC's VncViewer. This is a java application, so you need to have java installed on your computer.
https://java.com/en/download/help/download_options.xml

Download TigerVNC's VncViewer-1.9.0.jar (or whichever is the latest release) to use as a client on your local machine from:
https://github.com/TigerVNC/tigervnc/releases

Documentation for VncViewer is at:
https://tigervnc.org/doc/vncviewer.html

Info on TigerVNC:
https://tigervnc.org/

----------
2. Set up an ssh tunnel to connect from your computer to the VNC server.

2a. From Windows:

On Windows, you need to download an ssh client. PuTTY is the most common free client:
https://www.chiark.greenend.org.uk/~sgtatham/putty/

After downloading the client, open it by double-clicking the .exe executable file, then:
- In the tree on the left, navigate to Connection->SSH->Tunnels. 
    - Under "Add new forwarded port":
        - For the "Source Port", enter 5901 
        - For the "Destination Port", enter the floating IP address that you 
          assigned to your OpenStack server when you created it, followed by :5901 . If you can't 
          remember the floating ip address, then from within your OpenStack session, enter the command:
              openstack server list
          and the floating ip will be the second ip listed under Networks for your server. It 
          should look something like 10.20.0.XXX with the XXX replaced by some specific number. So 
          the final entry in your Destionation Port field in PuTTY should look something like:
              10.20.0.XXX:5901
        - Make sure the "Local" radio button is selected, and either "Auto" or "IPv4".
        - Click Add.
- Now, in the tree on the left, navigate to the Session section in PuTTY.
    - Under "Host Name (or IP address)", enter this, but replace "username" with your borealcloud username:
        username@borealcloud.nfis.org
    - Under "Port", enter 22
    - Make sure "SSH" is selected in "Connection type".
    - Click "Save".
    - Click "Open" to open an ssh session to the server. The password you will need to enter is your
      borealcloud account password.

2b. From Linux:

From a computer within PFC:
From a Terminal, do either of these two commands (I prefer the second one if I don't want to 
run the ssh tunneling in the background), but first replace "username" with your borealcloud 
username, and 10.20.0.XXX with the floating ip you assigned to your VNC server when you 
created it (see *below if you forget what that ip is):
    ssh -v -L 5901:10.20.0.XXX:5901 -N -f -l username borealcloud.nfis.org 
    ssh -v -L 5901:10.20.0.XXX:5901 -N -l username borealcloud.nfis.org 
The password you will need to enter is your borealcloud account password.

To connect from your home Linux computer, you have to use two tunnels: 
    ssh -v -L 2223:borealjump.nfis.org:22 -N -l brancourt borealremote.nfis.org
Then do a second ssh tunnel through the first one for VNC using either:
    ssh -v -L 5901:10.20.0.120:5901 -N -p 2223 brancourt@localhost
    ssh -v -L 5901:10.20.0.120:5901 -N -p 2223 -l brancourt localhost

*If you have forgotten what the floating ip address is for your VNC server, enter this command 
from within an OpenStack session:
    openstack server list
and the floating ip will be the second ip listed under Networks for your server. It should look 
something like 10.20.0.XXX with the XXX replaced by some specific number.

----------
3. Connect to your VNC server with TigerVNC's VncViewer:

3a. From Windows:

Open a Windows Command Prompt by clicking Start and then typing "cmd" into the search box.
Navigate to the folder that has the VncViewer-1.9.0.jar file which you downloaded in step 1. 
Then, launch the VncViewer client with this command:
    java -jar VncViewer-1.9.0.jar

In the connection details window that pops up, under "VNC server:" enter: localhost:5901
Then click Connect.
In the Password window that pops up, enter: centos
You should now be logged into your VNC server, showing graphical window displaying the Desktop.


3b. From Linux:

Open a new Terminal window and navigate to the folder that has the VncViewer-1.9.0.jar file 
which you downloaded in step 1. Then, launch the VncViewer client with this command:
    java -jar VncViewer-1.11.0.jar

In the connection details window that pops up, under "VNC server:" enter: localhost:5901
Then click Connect.
In the Password window that pops up, enter: centos
You should now be logged into your VNC server, showing graphical window displaying the Desktop.


----------
VNC Viewer Troubleshooting stuff...


If you get this error when trying to use VncViewer:
com.tigervnc.rdr.SystemException: javax.net.ssl.SSLHandshakeException: No appropriate protocol (protocol is disabled or cipher suites are inappropriate)

I think I solved this by going to Options-->Security and deselecting "Anonymous TLS" and selecting "None".

----------

----------

================================================================
Info on connecting via ssh tunneling:
https://unix.stackexchange.com/questions/115897/whats-ssh-port-forwarding-and-whats-the-difference-between-ssh-local-and-remot
http://www.allgoodbits.org/articles/view/58
https://community.time4vps.com/discussion/120/how-to-connect-to-vnc-server-ussing-ssh-tunnel
http://www.mit.edu/~avp/lqcd/ssh-vnc.html

--
To connect from a computer in the PFC building:
To connect with vncserver via the TigerVNC VncViewer java client, on your client machine, first do either of (I prefer the second if I don't want to run the ssh tunneling in the background):
ssh -v -C -L 5901:10.20.0.120:5901 -N -f -l brancourt borealcloud.nfis.org 
ssh -v -C -L 5901:10.20.0.120:5901 -N -l brancourt borealcloud.nfis.org 

Important! Make sure the ip address above is the ip of the OpenStack server that you connect to. 
The password that you will be prompted to enter at this stage is the password of your account on borealcloud.nfis.org. 

Then, launch the VncViewer client from a new terminal:
java -jar VncViewer-11.0.jar

Then connect to localhost:5901 (which connects to your local end of the ssh tunnel you set up, and forwards it on).
The password that you will be prompted to enter is the passord I set up for the vncserver, which is just "centos", the same as the username on the server.

--
To connect from your home Linux computer, you have to use two tunnels: 
ssh -v -L 2223:borealjump.nfis.org:22 -N -l brancourt borealremote.nfis.org
Then do a second ssh tunnel through the first one for VNC using either:
ssh -v -L 5901:10.20.0.120:5901 -N -p 2223 brancourt@localhost
ssh -v -L 5901:10.20.0.120:5901 -N -p 2223 -l brancourt localhost

Then, launch the VncViewer client from a new terminal:
java -jar VncViewer-11.0.jar

Then connect to localhost:5901 (which connects to your local end of the ssh tunnel you set up, and forwards it on).
The password that you will be prompted to enter is the passord I set up for the vncserver, which is just "centos", the same as the username on the server.

--
-L Specifies that the given port on the local (client) host is to be forwarded to the given host and port on the remote side.
ssh -L sourcePort:forwardToHost:onPort connectToHost means: connect with ssh to connectToHost, and forward all connection attempts to the local sourcePort to port onPort on the machine called forwardToHost, which can be reached from the connectToHost machine.

-C Requests compression of all data (including stdin, stdout, stderr, and data for forwarded X11, TCP and UNIX-domain connections).  The compression algorithm is the same used by gzip(1).  Compression is desirable on modem lines and other slow connections, but will only slow down things on fast networks. 
-N specifies that we're merely forwarding traffic not running a command on the remote machine.
-T disables pseudo-tty allocation.
-L indicates that the first port is a local port
-f : Requests ssh to go to background just before command execution. Once password supplied it will go to background and you can use prompt for type commands on local system.
-l vncuser : vncuser is the user to log in as on the remote machine.


==========
VNC Troubleshooting stuff...


Troubleshooting the error that occured when trying to use VncViewer on my home CentOS 8 machine:
com.tigervnc.rdr.SystemException: javax.net.ssl.SSLHandshakeException: No appropriate protocol (protocol is disabled or cipher suites are inappropriate)

CConn: connected to host localhost port 5901
CConnection: Server supports RFB protocol version 3.8
CConnection: Using RFB protocol version 3.8
com.tigervnc.rdr.SystemException: javax.net.ssl.SSLHandshakeException: No appropriate protocol (protocol is disabled or cipher suites are inappropriate)
	at com.tigervnc.rfb.CSecurityTLS.processMsg(CSecurityTLS.java:145)
	at com.tigervnc.rfb.CSecurityStack.processMsg(CSecurityStack.java:41)
	at com.tigervnc.rfb.CSecurityVeNCrypt.processMsg(CSecurityVeNCrypt.java:177)
	at com.tigervnc.rfb.CConnection.processSecurityMsg(CConnection.java:263)
	at com.tigervnc.rfb.CConnection.processMsg(CConnection.java:133)
	at com.tigervnc.vncviewer.VncViewer.run(VncViewer.java:565)
	at java.base/java.lang.Thread.run(Thread.java:834)
[ben@localhost TigerVNC]$ 
[ben@localhost TigerVNC]$ 

I think I solved this by going to Options-->Security and deselecting "Anonymous TLS" and selecting "None".


