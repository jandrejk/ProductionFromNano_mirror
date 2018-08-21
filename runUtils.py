import glob
import subprocess as sp
import shutil
import shlex
import json
import os

def main():
    # checkTokens()
    # checkProxy()
    getHeplxPublicFolder()

def checkTokens():

        neededToken = getSystem(inverse=True)
        if neededToken == "cern.ch": return True

        p = sp.Popen( shlex.split("klist"), stdout=sp.PIPE, stderr=sp.PIPE )
        (out,err) =  p.communicate()

        token = ""
        for line in out.splitlines():
            if "Ticket cache" in line:
                token = line.split("FILE:")[1]
            if "Default principal" in line:
                if not neededToken in line.lower():
                    print "Get a kerberos token for: {0}".format( neededToken )
                    return False
                else:
                    if not os.path.exists("kerberos"):
                        os.mkdir("kerberos")
                    shutil.copyfile(token, "kerberos/krb5_token")
                    return True

def checkProxy():

        p = sp.Popen( shlex.split("voms-proxy-info"), stdout=sp.PIPE, stderr=sp.PIPE )
        (out,err) =  p.communicate()

        if err:
            print err
            return False

        info = {"time":-1,"path":""}
        for line in out.splitlines():
            if "timeleft" in line: info['time'] = int( line.replace(" ","").split("timeleft:")[1].replace(":","") )
            if "path" in line: info["path"] = line.split(":")[1].replace(" ","")

        if info['time'] > 0:
            if not os.path.exists("proxy"):
                os.mkdir("proxy")
            shutil.copyfile(info["path"], "proxy/x509_proxy")
            return True
        else:
            print "Proxy expired! Get a new one with 'voms-proxy-init --voms cms'"
            return False

def getHeplxPublicFolder():

    user = os.environ["USER"]
    first = user[0]
    second = user[:8]

    return glob.glob( "/afs/hephy.at//user/{0}/{1}*/public".format(first,second) )[0]

def getSystem(inverse = False):
    host = os.environ["HOSTNAME"]
    if "heplx" in host: 
      if inverse: return "cern.ch"
      else: return "hephy.at"

    elif "lxplus" in host: 
      if inverse: return "hephy.at"
      else: return "cern.ch"

if __name__ == '__main__':
    main()