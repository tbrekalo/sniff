.PHONY: \
	all \
	debug \
	relwithdebinfo \
	release \
	clean

GEN:="Unix Makefiles"
NINJA_PATH:=$(shell command -v ninja 2> /dev/null)
ifdef NINJA_PATH
	GEN:="Ninja"
endif

build-debug: conanfile.txt
	conan install . --output-folder=$@ \
		--build=missing \
		-s build_type=Debug -s compiler.cppstd=gnu20; \
	cd $@; \
	cmake .. -G $(GEN) -DCMAKE_BUILD_TYPE=Debug \
		-DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake;

debug: build-debug
	cd build-debug; \
	cmake --build .;

clean-debug:
	rm -rf build-debug;

build-relwithdebinfo: conanfile.txt
	conan install . --output-folder=$@ \
		--build=missing \
		-s build_type=RelWithDebInfo -s compiler.cppstd=gnu20; \
	cd $@; \
	cmake .. -G $(GEN) -DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake;

relwithdebinfo: build-relwithdebinfo
	cd build-relwithdebinfo; \
	cmake --build .;

clean-relwithdebinfo:
	rm -rf build-relwithdebinfo;

build: conanfile.txt
	conan install . --output-folder=$@ \
		--build=missing \
		-s build_type=Release -s compiler.cppstd=gnu20; \
	cmake -B $@ -G $(GEN) -DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake;

release: build
	cd build; \
	cmake --build .;

clean-release:
	rm -rf build;

clean: clean-debug clean-relwithdebinfo clean-release
	@:

all: debug relwithdebinfo release
	@:
