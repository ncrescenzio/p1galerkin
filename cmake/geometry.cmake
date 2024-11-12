# Search "geometry" library
if(NOT TARGET "geometry")

   # Try to find "geometry" library using
   # -Dgeometry_DIR, -Dgeometry_ROOT, -DGEOMETRY_DIR and -DGEOMETRY_ROOT
   find_package(geometry PATHS ${GEOMETRY_ROOT} ${GEOMETRY_DIR} QUIET)

   if(NOT geometry_FOUND)

      # Assume that source code for geometry library
      # is in the same directory of linalg library
      set(GEOMETRY_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../geometry)

      if(EXISTS ${GEOMETRY_SRC})

         # Assume the code is already compiled
         if(EXISTS ${GEOMETRY_SRC}/build)
            find_package(geometry REQUIRED
               PATHS ${GEOMETRY_SRC}/build/
               NO_DEFAULT_PATH)

         # The code is not compiled
         else(EXISTS ${GEOMETRY_SRC}/build)
            add_subdirectory(${GEOMETRY_SRC} ${GEOMETRY_SRC}/build EXCLUDE_FROM_ALL)

         endif(EXISTS ${GEOMETRY_SRC}/build)

      else(EXISTS ${GEOMETRY_SRC})
         message(FATAL_ERROR "Library 'geometry' NOT FOUND!")

      endif(EXISTS ${GEOMETRY_SRC})

   endif(NOT geometry_FOUND)

endif(NOT TARGET "geometry")
