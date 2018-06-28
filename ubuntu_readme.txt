Download Kakadu for Linux http://kakadusoftware.com/downloads/ and see Kakadu's README.txt for instructions with LD_LIBRARY_PATH. If you encounter the error "kakadu/kdu_compress: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.21' not found" do the following,

sudo add-apt-repository ppa:ubuntu-toolchain-r/test 
sudo apt-get update
sudo apt-get install 
sudo apt-get install libstdc++6


