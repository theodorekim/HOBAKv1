include ../include_top.mk

# calls:
INCLUDES   = ${INCLUDES_COMMON} -I../../ext/catch2
CFLAGS     = ${CFLAGS_COMMON} ${INCLUDES} -std=c++17 -DGL_DISABLED
LDFLAGS    = ${LDFLAGS_COMMON} -std=c++17
EXECUTABLE = ../../bin/regressionTests

SOURCES    = regressionTests.cpp \
						 ../../ext/catch2/catch_amalgamated.cpp \
						 ../../src/util/COLLISION_UTIL.cpp \
						 ../../src/util/MATRIX_UTIL.cpp \
						 ../../src/util/MIN_HEAP.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/util/FILE_IO.cpp \
						 ../../src/Timestepper/QUASISTATIC.cpp \
						 ../../src/Timestepper/NEWMARK.cpp \
						 ../../src/Timestepper/BACKWARD_EULER_VELOCITY.cpp \
						 ../../src/Timestepper/BACKWARD_EULER_POSITION.cpp \
						 ../../src/Timestepper/TIMESTEPPER.cpp \
						 ../../src/Geometry/LINE_INTERSECT.cpp \
						 ../../src/Geometry/FIELD_3D.cpp \
						 ../../src/Geometry/FIELD_2D.cpp \
						 ../../src/Geometry/CAPSULE.cpp \
						 ../../src/Geometry/CYLINDER.cpp \
						 ../../src/Geometry/CUBE.cpp \
						 ../../src/Geometry/SPHERE.cpp \
						 ../../src/Geometry/TET_MESH.cpp \
						 ../../src/Geometry/TET_MESH_FASTER.cpp \
						 ../../src/Geometry/AABB_TREE.cpp \
						 ../../src/Damping/Volume/GREEN_DAMPING.cpp \
						 ../../src/Hyperelastic/Volume/MCADAMS_COLLISION.cpp \
						 ../../src/Hyperelastic/Volume/EDGE_COLLISION.cpp \
						 ../../src/Hyperelastic/Volume/EDGE_SQRT_COLLISION.cpp \
						 ../../src/Hyperelastic/Volume/EDGE_HYBRID_COLLISION.cpp \
						 ../../src/Hyperelastic/Volume/VERTEX_FACE_COLLISION.cpp \
						 ../../src/Hyperelastic/Volume/VERTEX_FACE_SQRT_COLLISION.cpp \
						 ../../src/Hyperelastic/Volume/ANISOTROPIC_ARAP.cpp \
						 ../../src/Hyperelastic/Volume/ANISOTROPIC_STVK.cpp \
						 ../../src/Hyperelastic/Volume/ANISOTROPIC_FUNG.cpp \
						 ../../src/Hyperelastic/Volume/ANISOTROPIC_DIRICHLET.cpp \
						 ../../src/Hyperelastic/Volume/STVK.cpp \
						 ../../src/Hyperelastic/Volume/SNH_WITH_BARRIER.cpp \
						 ../../src/Hyperelastic/Volume/SNH.cpp \
						 ../../src/Hyperelastic/Volume/NEO_HOOKEAN_BW.cpp \
						 ../../src/Hyperelastic/Volume/ARAP.cpp \
						 ../../src/Hyperelastic/Volume/LINEAR.cpp \
						 ../../src/Hyperelastic/Volume/HYPERELASTIC.cpp

include ../include_bottom.mk
