import os,sys
def reject_original(event,args):
    if event=="open" and args and isinstance(args[0],(str,bytes,os.PathLike)):
        path=os.path.abspath(os.fsdecode(args[0]))
        if path=="/home/rihen/Desktop/ESR1" or path.startswith("/home/rihen/Desktop/ESR1/"):
            raise PermissionError("Standalone verification forbids original-repository access: "+path)
    if event=="socket.connect":
        raise PermissionError("Offline verification forbids network connections")
sys.addaudithook(reject_original)
