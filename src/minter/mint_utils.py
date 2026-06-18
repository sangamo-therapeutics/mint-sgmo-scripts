import os

try:
    HOST_UID = int(os.getenv("HOST_UID"))
except (TypeError, ValueError):
    HOST_UID = None
try:
    HOST_GID = int(os.getenv("HOST_GID"))
except (TypeError, ValueError):
    HOST_GID = None


def set_newfile_permissions(directory, prior_files=None, host_uid=HOST_UID, host_gid=HOST_GID, **kwargs):
    if not host_uid:
        return
    prior_files = prior_files or []
    new_files = directory.glob(f"*.*")
    new_files = [f for f in new_files if f not in prior_files]
    for nf in new_files:
        try:
            os.chown(nf, host_uid, host_gid)

            os.chown(nf, HOST_UID, HOST_GID)
        except Exception as e:
            print(f"unable to change permissions on file: {nf} {e}")
