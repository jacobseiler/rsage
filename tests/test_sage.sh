currdir=`pwd`
cd tests/

python3 test_sage.py
exit_code=$? 
if [ $exit_code -ne 0 ]; then
  echo "test_sage.py exited with errorcode $exit_code"
  exit $exit_code 
fi

cd "$currdir"
