/bin/docker build . --build-arg UID=$(id -u) -t parafem
cwd=$(pwd)
dir=$(dirname ${cwd})
echo ${dir}
docker run -t -v ${dir}:/home/conda/parafem parafem conda-build parafem --output-folder parafem/packages --python=3.7 --dirty