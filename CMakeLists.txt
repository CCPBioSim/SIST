cmake_minimum_required(VERSION 3.28)
project(code)

set(CMAKE_CXX_STANDARD 14)

include_directories(src/trans_compete)
include_directories(src/trans_three)

add_executable(code
        src/trans_compete/G_x.cpp
        src/trans_compete/G_x.h
        src/trans_compete/qsidd.cpp
        src/trans_compete/SIDD_1R.cpp
        src/trans_compete/SIDD_1R.h
        src/trans_compete/SIDD_2R.cpp
        src/trans_compete/SIDD_2R.h
        src/trans_compete/SIDD_3R.cpp
        src/trans_compete/SIDD_3R.h
        src/trans_compete/SIDD_4R.cpp
        src/trans_compete/SIDD_4R.h
        src/trans_compete/SIDD_Base.cpp
        src/trans_compete/SIDD_Base.h
        src/trans_compete/stat_1R.cpp
        src/trans_compete/stat_1R.h
        src/trans_three/G_x.cpp
        src/trans_three/G_x.h
        src/trans_three/qsidd.cpp
        src/trans_three/SIDD_1R.cpp
        src/trans_three/SIDD_1R.h
        src/trans_three/SIDD_2R.cpp
        src/trans_three/SIDD_2R.h
        src/trans_three/SIDD_3R.cpp
        src/trans_three/SIDD_3R.h
        src/trans_three/SIDD_4R.cpp
        src/trans_three/SIDD_4R.h
        src/trans_three/SIDD_Base.cpp
        src/trans_three/SIDD_Base.h
        src/trans_three/stat_1R.cpp
        src/trans_three/stat_1R.h)
