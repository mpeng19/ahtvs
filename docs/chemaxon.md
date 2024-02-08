
ChemAxon Setup OSX (THIS IS FOR MEMBERS OF ASPURU-GUZIK GROUP ONLY)
-------------------

10. head to [ChemAxon](https://www.chemaxon.com) and register on their website.  This is the normal name, password, email, activate process.  You can make your own login, there is no shared login.

20. Download [JChem Suite for OS X](https://www.chemaxon.com/download.php?d=/data/download/jchem/16.6.20.0/jchem-16.6.20.0-macos.dmg&dl=true)

30. Run that installer

40. Run `open -a LicenseManager.app`

50. Browse to `a2g2/docs/license_2014.cxl`  and click install.



ChemAxon Setup LINUX (THIS IS FOR MEMBERS OF ASPURU-GUZIK GROUP ONLY)
---------------------

10. head to [ChemAxon](https://www.chemaxon.com) and register on their website.  This is the normal name, password, email, activate process.  You can make your own login, there is no shared login.

20. Download [JChem Suite for Linux 64 bit](https://www.chemaxon.com/download.php?d=/data/download/jchem/6.1.5/jchem-6.1.5-linux_with_jre64.sh)

30. execute the above downloaded shell script (use as sudo and install the software into /opt/ using their gui installer):

        chmod 755 jchem-6.1.5-linux_with_jre64.sh
        sudo ./jchem-6.1.5-linux_with_jre64.sh

40. install java (lifted from http://ubuntuhandbook.org/index.php/2013/07/install-oracle-java-6-7-8-on-ubuntu-13-10/)

        sudo add-apt-repository ppa:webupd8team/java
        sudo apt-get update
        sudo apt-get install oracle-java8-installer

50. add the ChemAxon bin dirs to your path:  append these lines to your .bashrc

        export PATH=$PATH:/opt/ChemAxon/JChem/bin/
        export PATH=$PATH:/opt/ChemAxon/MarvinBeans/bin/

60. now open new terminal or run `source ~/.bashrc`

65. install license:

    a. run `/opt/ChemAxon/JChem/bin/license`
    b. click browse and select file in this repo: `a2g2/docs/license_2014.cxl`
    c. click install button

70. for more help use -h on apps:

        mview -h
        react -h