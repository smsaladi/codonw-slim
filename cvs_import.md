Importing from CVS
==================

git cvsimport -v -a -k -A authors.txt -d `pwd`/cvs -C codonw Src
git cvsimport -v -a -k -A authors.txt -d `pwd`/cvs -C codonw_docs Documentation


cd codonw
git branch -D origin
git branch -m cvsimport
git checkout -b master

git remote add docs ../codonw_docs
git checkout -b gh-pages docs/master

git remote rm docs
rm -rf ../codonw_docs

