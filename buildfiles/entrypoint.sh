#!/bin/bash -il

# use non-interactive cp
unalias cp
# Create conda user with the same uid as the host, so the container can write
# to mounted volumes
# Adapted from https://denibertovic.com/posts/handling-permissions-with-docker-volumes/
USER_ID=${HOST_USER_ID:-9001}
echo "User_id: ${USER_ID}"
useradd --shell /bin/bash -u $USER_ID -G lucky -o -c "" -m conda
export HOME=/home/conda
export USER=conda
export LOGNAME=conda
export MAIL=/var/spool/mail/conda
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/conda/bin
if [ "$(uname -m)" = "x86_64" ]; then
   export supkg="gosu"
elif [ "$(uname -m)" = "aarch64" ]; then
   export supkg="su-exec"
else
   export supkg="su-exec"
fi
chown conda:conda $HOME
cp -Rn /etc/skel $HOME && chown -R conda:conda $HOME/skel && (ls -A1 $HOME/skel | xargs -I {} mv -fn $HOME/skel/{} $HOME) && rm -Rf $HOME/skel
cp -n /root/.condarc $HOME/.condarc && chown conda:conda $HOME/.condarc
cd $HOME

# Source everything that needs to be.
. /opt/docker/bin/entrypoint_source

# Run whatever the user wants.
exec /opt/conda/bin/$supkg conda "$@"