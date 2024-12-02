cmake_minimum_required(VERSION 3.0)

set(CXX_STANDARD_REQUIRED On)
set(CXX_STANDARD C++11)

set(LIBRARY_PATH_EIGEN "../../../libraries/eigen/")

# Spécifie les options passées au compilateur. Ici : l'optimisation du code binaire et les warnings.
set(CMAKE_CXX_FLAGS "-Wall")
# set(CMAKE_CXX_FLAGS_DEBUG "-g") # Les paramètres par défaut: décommentez si vous voulez les changer
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG") # Les paramètres par défaut: décommentez si vous voulez les changer

# Options diverses de CMake, à ne pas modifier.
set(CMAKE_SKIP_INSTALL_RULES True)

# Décommentez la ligne que vous souhaitez, puis recompilez.
set(CMAKE_BUILD_TYPE "Debug")
#set(CMAKE_BUILD_TYPE "Release")

# Options spécifiques à ce projet.
project(CHP)

# Liste exhaustive des fichiers à compiler ensemble.
set(PROJECT_COMPILATION_FILE_LIST read.cpp fonctions.cpp Matrice.cpp BiCGstab.cpp)

# Précise que l'exécutable sera à assembler avec ces fichiers compilés.
add_executable(${PROJECT_NAME} ${PROJECT_COMPILATION_FILE_LIST})

# Ajoute les répertoires des bibliothèques liées, ici Eigen et les petites biblio locales au projet.
target_include_directories(${PROJECT_NAME} PUBLIC "${LIBRARY_PATH_EIGEN}" "src/libraries")

# Spécifie les bibliothèques MPI à inclure lors de la liaison du programme.
find_package(MPI REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)
