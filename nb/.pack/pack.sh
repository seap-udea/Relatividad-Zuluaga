#!/bin/bash
basedir=$(pwd)
pack=$1
if [ "x$pack" == "x" ];then pack="pack";fi

storedir=".pack"
confile="pack.conf"

if [ $pack = "pack" ];then
    echo "Packing..."
    rm -r $storedir/*--* &> /dev/null
    for file in $(cat $storedir/$confile |grep -v "#")
    do
    	echo -e "\tfile $file..."
    	fname=$(basename $file)
    	dname=$(dirname $file)
    	uname=$(echo $dname |sed -e s/\\//_/)
    	sdir="$storedir/$uname--$fname"
    	mkdir -p "$sdir"
    	cd $sdir
    	split -b 20000k $basedir/$file $fname-
    	cd - &> /dev/null
	if [ -d .git ];then 
	    echo "$fname" >> $dname/.gitignore
    	    git add "$storedir/$uname--$fname/*"
	fi
    done
    if [ -d .git ];then 
	git add $storedir/*--*
    fi
else
    echo "Unpacking..."
    for file in $(cat $storedir/$confile |grep -v "#")
    do
    	echo -en "\tUnpacking $file..."
    	fname=$(basename $file)
    	dname=$(dirname $file)
    	uname=$(echo $dname |sed -e s/\\//_/)
    	sdir="$storedir/$uname--$fname"
    	if [ -e $dname/$fname ];then
    	    echo "file already unpacked"
    	else
    	    cat "$sdir"/$fname-* > $dname/$fname
	    echo "$fname" >> $dname/.gitignore	    
	    if [ -d .git ];then 
    		git add "$dname/.gitignore"
	    fi
    	    echo
    	fi
    done
fi
