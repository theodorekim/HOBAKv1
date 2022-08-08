include ../include_top.mk

# calls:
INCLUDES   = ${INCLUDES_COMMON} -I../../ext/glvu/ -I../../src/Geometry/
CFLAGS     = ${CFLAGS_COMMON} ${INCLUDES}
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = ../../bin/closestPointViewer3D

SOURCES    = closestPointViewer3D.cpp \
						../../src/Geometry/FIELD_2D.cpp \
						../../src/util/TIMER.cpp \
						../../src/util/MIN_HEAP.cpp \
						../../src/Geometry/VECTOR3_FIELD_3D.cpp \
						../../src/Geometry/FIELD_3D.cpp

include ../include_bottom.mk
