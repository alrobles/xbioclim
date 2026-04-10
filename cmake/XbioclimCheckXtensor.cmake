# XbioclimCheckXtensor.cmake — verify xtensor is available with a helpful message
#
# Include this module after find_package(xtensor ...) to get a clear error
# message pointing users at the right install procedure.

if(NOT xtensor_FOUND)
    message(FATAL_ERROR
        "xtensor (>= 0.25) not found.\n"
        "Install from source:\n"
        "  1. Install xtl first:\n"
        "     git clone --depth=1 --branch 0.7.7 https://github.com/xtensor-stack/xtl.git /tmp/xtl\n"
        "     cmake -S /tmp/xtl -B /tmp/xtl/build -DCMAKE_INSTALL_PREFIX=/usr/local\n"
        "     sudo cmake --build /tmp/xtl/build --target install\n"
        "  2. Then xtensor:\n"
        "     git clone --depth=1 --branch 0.25.0 https://github.com/xtensor-stack/xtensor.git /tmp/xtensor\n"
        "     cmake -S /tmp/xtensor -B /tmp/xtensor/build -DCMAKE_INSTALL_PREFIX=/usr/local\n"
        "     sudo cmake --build /tmp/xtensor/build --target install\n"
        "Or set xtensor_DIR / CMAKE_PREFIX_PATH to a custom installation.")
endif()
