# imp-server

This archival repository does not contain the current Beagle software.  The
current Beagle release can be obtained from
[Beagle 5.5 web page](https://faculty.washington.edu/browning/beagle/beagle.html).

You can create the Beagle jar file with the commands:

    git clone https://github.com/tmp-sharing/imp-server.git
    javac -cp imp-server/src/ imp-server/src/main/Main.java
    jar cfe beagle.jar main/Main -C imp-server/src/ ./

[Beagle 5.5 documentation](https://faculty.washington.edu/browning/beagle/beagle_5.5_17Dec24.pdf)

