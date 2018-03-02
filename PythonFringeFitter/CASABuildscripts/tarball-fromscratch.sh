#!/bin/bash
#-----------------------------------------------------------------------
#-- create a stand-alone tarball from a local casa build
#--
#-- version 20170607, by drr
#-- version 20170608, by drr
#-----------------------------------------------------------------------

#-- casa build results
srcdir=~/casa/pkg
#-- target directory
tgdir=~/casa/release/casa-dev-20170511-3
#-- official casa release
refdir=~/casa/release/casa-release-4.7.2-el7

#============
#-- packages
#============
rm -rf $tgdir
mkdir $tgdir
pushd $tgdir
dirs=(casacore code gcwrap asap)
for dir in "${dirs[@]}"; do
	pushd $srcdir/$dir || exit 1
	tar -cf - * | tar -xf - -C $tgdir
	popd
done
#-- move stuff around
mv python/2.7/* lib/python2.7/
rmdir python/2.7 python
mkdir -p lib/casa/bin
mv bin/* lib/casa/bin/
mv lib/casa/bin/qcasabrowser lib/casa/bin/casabrowser
cp -a $srcdir/etc .
rm -rf include
rm -rf share
popd


#============
#-- augment from casa RPMs, system, or official release
#============
#-- bin
cp -v /opt/casa/02/bin/{i,}python $tgdir/lib/casa/bin/
#cp -iva $refdir/plugins $tgdir/  #-- needed?
cp -iv --no-dereference $refdir/bin/* $tgdir/bin/
pushd $tgdir/bin
ln -sf ../lib/casa/bin/{,i}python .
popd
perl -ni -e 'print; print "} else {\n    \$ENV{CASALD_LIBRARY_PATH} = \"\$installpath/lib\";\n    \$ENV{LD_LIBRARY_PATH} = \"\$installpath/lib\";\n" if /\$ENV\{LD_LIBRARY_PATH\} = \$ENV\{CASALD_LIBRARY_PATH\}/' $tgdir/bin/casa
perl -pi -e 's+casa3party="/usr/lib/casapy"+casa3party=\${CASAPATH% *}+; s/installed_from_rpm="F"/installed_from_rpm="T"/' $tgdir/lib/casa/bin/casa
sed -i -e '50,53d' $tgdir/lib/casa/bin/casa

#-- python
cp -ia /opt/casa/02/lib/python2.7/* $tgdir/lib/python2.7/

#-- share
cp -ia $refdir/share $tgdir/

#-- lib
#-- collect file names
files=()
for f in $(\find $refdir -not \( -path "$refdir/etc/carta" -prune \) -name lib\*.so\*); do
	file=${f##*/}
	file=${file%%.so*}
	files+=("$file")
done
files+=("libselinux")

#-- uniq
files=($(echo "${files[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

#-- copy files
for file in "${files[@]}"; do
	ls $tgdir/lib/$file.so* &>/dev/null && continue
	echo
	echo $file.so
	if ls /opt/casa/02/lib/$file.so* &>/dev/null; then
		src=($(ls /opt/casa/02/lib/$file.so*))
	elif ls /usr/lib64/$file.so* &>/dev/null; then
		src=($(ls /usr/lib64/$file.so*))
	else
		src=($(locate "${file}.so" |sed -ne '/casa-release-4.7.2-el7/d; /NIFPGA/d; /casa-dev-/d; /casa\/01/d; p'))
	fi

	if [[ -n "${src[@]}" ]]; then
		cp -vi --no-dereference "${src[@]}" "$tgdir/lib/"
	fi
done


#============
#-- data
#============
excludefile=/tmp/exclude.txt
cat << EOF > $excludefile
regression/
.svn/
atnf/
protopipe/
demo/tutorial/
demo/benchmark/
demo/autoflag/
demo/jupiter6cm.fits
demo/atca4800.fits
demo/dishdemo/
demo/tutorial
alma/test/
alma/responses/ALMA_AIRY_12M.VP/
alma/responses/ALMA_AIRY_7M.VP/
nrao/VLA/*.fits
nrao/VLA/vlafiller_test
EOF
rsync -au /opt/casa/data --exclude-from=$excludefile $tgdir

rm $excludefile


#============
#-- tarball
#============
dirname=${tgdir##*/}
echo
pushd $tgdir/..
echo "tarring up $PWD/$tgdir.tar.xz ..."
tar -I pxz -cf $dirname.tar.xz $dirname
md5sum $dirname.tar.xz > $dirname.tar.xz.md5
rsync -va $dirname.tar.xz{,.md5} lilo:www/
popd

echo "done"
