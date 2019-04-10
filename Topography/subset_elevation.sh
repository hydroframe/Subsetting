#!/bin/sh
#require python3 and these libraries: gdal, ogr, osr, shapely, numpy
# yum install wget
# yum install openssl-devel
# yum install libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev
# yum -y groupinstall development
# yum -y install zlib-devel
# wget https://www.python.org/ftp/python/3.6.3/Python-3.6.3.tar.xz
# tar xJf Python-3.6.3.tar.xz
# cd Python-3.6.3
# ./configure
# make
# make install

###Parse inputs
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -sel_type)
    sel_type="$2"
    shift # past argument
    shift # past value
    ;;
    -s)
    if [ "$sel_type" = '-id' ]; then
    	s="$2"
    else
    	s="$2 $3"
    fi
    shift # past argument
    shift # past value
    ;;
    -out_name)
    out_name="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}"



###Download and install pf-mask-utilities
if [ ! -d pf-mask-utilities ]; then
	git clone https://github.com/smithsg84/pf-mask-utilities.git
fi

cd pf-mask-utilities
make

cd ..

chmod 755 subset_domain_by_shape.py
echo ./subset_domain_by_shape.py $sel_type $s $out_name

###Run python file
./subset_domain_by_shape.py $sel_type $s $out_name
