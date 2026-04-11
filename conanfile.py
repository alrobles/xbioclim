from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout


class XbioclimConan(ConanFile):
    name = "xbioclim"
    version = "0.1.0"
    license = "MIT"
    url = "https://github.com/alrobles/xbioclim"
    description = (
        "C++17 library for computing BIO01–BIO19 bioclimatic variables "
        "from monthly climate rasters using xtensor and GDAL."
    )
    topics = ("bioclimatic", "climate", "gdal", "xtensor", "geospatial")

    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "with_openmp": [True, False],
        "build_tests": [True, False],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "with_openmp": True,
        "build_tests": False,
    }

    exports_sources = (
        "CMakeLists.txt",
        "cmake/*",
        "include/*",
        "src/*",
        "tests/*",
        "tools/*",
    )

    def requirements(self):
        self.requires("gdal/3.6.2")
        self.requires("xtensor/0.24.7")
        if self.options.build_tests:
            self.requires("catch2/3.4.0")

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["XBIOCLIM_BUILD_TESTS"] = self.options.build_tests
        tc.variables["XBIOCLIM_USE_OPENMP"] = self.options.with_openmp
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["xbioclim"]
        if self.options.with_openmp:
            if self.settings.os in ("Linux", "FreeBSD"):
                self.cpp_info.system_libs.append("gomp")
