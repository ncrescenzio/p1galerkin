# Search "linalg" library
if(NOT TARGET "linalg")

   # Try to find "linalg" library using
   # -Dlinalg_DIR, -Dlinalg_ROOT, -DLINALG_DIR and -DLINALG_ROOT
   find_package(linalg PATHS ${LINALG_ROOT} ${LINALG_DIR} QUIET)

   if(NOT linalg_FOUND)

      # Assume that source code for linalg library
      # is in the same directory of linalg library
      set(LINALG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../linear_algebra)

      if(EXISTS ${LINALG_SRC})

         # Assume the code is already compiled
         if(EXISTS ${LINALG_SRC}/build)
            find_package(linalg REQUIRED
               PATHS ${LINALG_SRC}/build/
               NO_DEFAULT_PATH)

         # The code is not compiled
         else(EXISTS ${LINALG_SRC}/build)
            add_subdirectory(${LINALG_SRC} ${LINALG_SRC}/build EXCLUDE_FROM_ALL)

         endif(EXISTS ${LINALG_SRC}/build)

      else(EXISTS ${LINALG_SRC})
         message(FATAL_ERROR "Library 'linalg' NOT FOUND!")

      endif(EXISTS ${LINALG_SRC})

   endif(NOT linalg_FOUND)

endif(NOT TARGET "linalg")
