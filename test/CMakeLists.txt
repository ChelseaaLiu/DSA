cmake_minimum_required(VERSION 3.20)

project(TestDSA CXX)

file(GLOB files CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR} *.cpp *.h)

# message(STATUS "files: ${files}")
add_executable(${PROJECT_NAME} ${files})
target_include_directories(${PROJECT_NAME} PUBLIC .)
target_link_libraries(${PROJECT_NAME}
    PUBLIC
    DSA
    PUBLIC
    gtest_main
)
target_compile_definitions(${PROJECT_NAME} PRIVATE PROJECTDIR="${PROJECT_SOURCE_DIR}")