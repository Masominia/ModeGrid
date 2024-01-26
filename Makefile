.PHONY: install build test lint run all

install:
	pip install -e ".[dev]"

CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Icpp/include
LIB := $(filter-out cpp/src/main.cpp,$(wildcard cpp/src/*.cpp))

build:
	mkdir -p cpp/build
	$(CXX) $(CXXFLAGS) -o cpp/build/modegrid cpp/src/main.cpp $(LIB)
	$(CXX) $(CXXFLAGS) -o cpp/build/modegrid_tests $(wildcard cpp/tests/*.cpp) $(LIB)

lint:
	ruff check src tests scripts
	ruff format --check src tests scripts

test: build
	cpp/build/modegrid_tests
	pytest

run: build
	python -m modegrid run

all: lint test run
