export PKG_CONFIG_PATH=/home/shuaiw/miniconda3/envs/cpp/lib/pkgconfig/
export LD_LIBRARY_PATH=/home/shuaiw/miniconda3/envs/cpp/lib:$LD_LIBRARY_PATH
g++ -o test load_bam.cpp $(pkg-config --cflags --libs htslib pbbam)
