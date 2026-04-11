from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout
from conan.tools.files import copy
import os


class XbioclimConan(ConanFile):
    name = "xbioclim"
    version = "0.1.0"
    description = (
        "C++17 library for computing BIO01-BIO19 bioclimatic variables "
        "from monthly climate rasters using GDAL and xtensor"
    )
    license = "MIT"
    url = "https://github.com/alrobles/xbioclim"
    homepage = "https://github.com/alrobles/xbioclim"
    topics = ("bioclimate", "gdal", "xtensor", "raster", "climate")

    settings = "os", "compiler", "build_type", "arch"

    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "with_openmp": [True, False],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "with_openmp": True,
    }

    exports_sources = (
        "CMakeLists.txt",
        "cmake/*",
        "include/*",
        "src/*",
    )

    def requirements(self):
        self.requires("gdal/[>=3.0 <4]")
        self.requires("xtensor/0.25.0")

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["XBIOCLIM_BUILD_TESTS"] = False
        tc.variables["XBIOCLIM_USE_OPENMP"] = self.options.with_openmp
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(self, "LICENSE", self.source_folder, os.path.join(self.package_folder, "licenses"))
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["xbioclim"]
        self.cpp_info.set_property("cmake_target_name", "xbioclim::xbioclim")
        if self.options.with_openmp:
            if self.settings.os == "Linux":
                self.cpp_info.system_libs.append("gomp")
            elif self.settings.compiler == "apple-clang":
                self.cpp_info.system_libs.append("omp")
