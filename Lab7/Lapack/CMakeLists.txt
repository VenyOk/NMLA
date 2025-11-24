cmake_minimum_required(VERSION 3.16)
project(Lapack)

# Поиск необходимых пакетов
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Включаемые директории для LAPACK
include_directories(${LAPACK_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} main.cpp)

# Линковка библиотек
target_link_libraries(${PROJECT_NAME} PRIVATE 
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

# Для Linux/Unix систем может понадобиться добавить -llapacke
if(UNIX)
    target_link_libraries(${PROJECT_NAME} PRIVATE lapacke)
endif()