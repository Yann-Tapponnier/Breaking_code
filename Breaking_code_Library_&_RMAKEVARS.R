###########################################################################################################################
######################################################## RMAKEVARS ########################################################
###########################################################################################################################


Use to direct the R sessions scripts to look at the right place for the compilators
##################################################################################

# In R :
usethis::edit_r_makevars()
#or 
nano ~/.R/Makevars

                        # #Old RMakevars Enter all of that : 
                        # LOC=/opt/homebrew/Cellar/gcc/15.2.0_1
                        # 
                        # CC=$(LOC)/bin/gcc-15 -fopenmp
                        # CXX=$(LOC)/bin/g++-15 -fopenmp
                        # CXX11=$(LOC)/bin/g++-15 -fopenmp
                        # 
                        # CFLAGS=-g -O3 -Wall -pedantic -std=gnu99 -mtune=generic -pipe
                        # CXXFLAGS=-g -O3 -Wall -pedantic -std=c++11 -mtune=generic -pipe
                        # # LDFLAGS=-L$(LOC)/lib -Wl,-rpath,$(LOC)/lib,-L/usr/local/lib
                        # # CPPFLAGS=-I$(LOC)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -I/usr/local/include -I/opt/homebrew/opt/boost/include
                        # 
                        # FLIBS=-L$(LOC) -L$(LOC)/lib -lgfortran -lquadmath -lm
                        # CXX1X=$(LOC)/bin/g++-15
                        # CXX98=$(LOC)/bin/g++-15
                        # CXX11=$(LOC)/bin/g++-15
                        # CXX14=$(LOC)/bin/g++-15
                        # CXX17=$(LOC)/bin/g++-15
                        # 
                        # 
                        # LDFLAGS=-L$(LOC)/lib -L/usr/local/lib -L/opt/homebrew/opt/zlib/lib -Wl,-rpath,$(LOC)/lib
                        # CPPFLAGS=-I$(LOC)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -I/usr/local/include -I/opt/homebrew/opt/boost/include -I/opt/homebrew/opt/zlib/include
                        # 

          
################## NEW RMAKEVARS
LOC=/opt/anaconda3/

CC=clang
CXX=clang++
LDFLAGS=-L$(LOC)/lib -L/usr/local/lib -L/opt/homebrew/opt/zlib/lib -Wl,-rpath,$(LOC)/lib
CPPFLAGS=-I$(LOC)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -I/usr/local/include -I/opt/homebrew/opt/boost/include -I/opt/homebrew/opt/z$
PKG_CONFIG_PATH="/opt/homebrew/opt/openssl/lib/pkgconfig:/opt/homebrew/opt/zlib/lib/pkgconfig"


F77=/opt/anaconda3/bin/gfortran
FC=/opt/anaconda3/bin/gfortran

CXXFLAGS=-O2 -Wall -std=c++11
CXX11FLAGS=-O2 -Wall -std=c++11
CXX14FLAGS=-O2 -Wall -std=c++14
CXX17FLAGS=-O2 -Wall -std=c++17




##### Rstudio does not rely on zshrc sor you need 
#launchctl setenv RSTUDIO_WHICH_R /opt/anaconda3/bin/R
launchctl setenv RSTUDIO_WHICH_R /usr/local/bin/R
#--> Copy it to Z/profile to be permanante :
nano ~/.zprofile


# in r-env
launchctl setenv RSTUDIO_WHICH_R /Users/administrateur/micromamba/envs/r_env/bin/R
          
#### You can revert this linking: 
launchctl unsetenv RSTUDIO_WHICH_R


          
# 
# ########## COPYING ALL THE PREVIOUS LIBRARIS
#           # It was all the libraries installed in a previous version of R
# /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
# 
# #Now I point to
# /opt/anaconda/R/and I copy paste the lib already intalled previously
# /opt/anaconda3/lib/R/library








