# C compiler and flags
CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-O3 -g

OBJECT = KHI.cpp 
HELPER = input.h
TARGET = KHI

all : $(TARGET)

$(TARGET) : $(OBJECT) $(HELPER)
	$(CXX) $(CPPFLAGS) -o $@ $(OBJECT)
