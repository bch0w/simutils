conda create -n pyatoa python=3.7 -y
source activate pyatoa
git clone --branch devel https://github.com/bch0w/pyatoa.git
cd pyatoa
pip install -e .
cd ..
git clone --branch devel https://github.com/bch0w/seisflows3.git
cd seisflows3
pip install -e .
