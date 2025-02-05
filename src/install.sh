export PKG_CONFIG_PATH=/home/shuaiw/miniconda3/envs/cpp/lib/pkgconfig/
g++ -o test load_bam.cpp $(pkg-config --cflags --libs htslib)
