SHELL := /bin/bash

OBJ_DIR = obj
SRC_DIR = ../../src
INC_DIR = ../../include
DAT_DIR = ../../data
LIB_DIR = ../

TST_DIR = ../../tests/*

LST_DIR = $(INC_DIR) \
		  $(DAT_DIR)

INCLUDE := $(foreach inc, $(LST_DIR), -I$(inc))

CC = g++
CCFLAGS = -g -ansi -std=c++17 $(INCLUDE)
LDFLAGS = -lmtqr -lm -lquadmath -lgsl -lgslcblas

OBJ = $(OBJ_DIR)/mtqr.o $(OBJ_DIR)/data_management.o $(OBJ_DIR)/monomial_transformation.o $(OBJ_DIR)/vector_operations.o

LIB = $(LIB_DIR)libmtqr.a