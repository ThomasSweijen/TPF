# File generated by kdevelop's qmake manager. 
# ------------------------------------------- 
# Subdir relative project main directory: ./plugins/Engine/Engine/InteractionEngine/NarrowInteractionGeometryEngine/SwiftPolyhedronProximityModeler
# Target is a library:  

LIBS += -lPolyhedralSweptSphere \
        -lMacroMicroContactGeometry \
        -lyade-lib-swiftpp \
        -rdynamic 
INCLUDEPATH += $(YADEINCLUDEPATH) 
MOC_DIR = $(YADECOMPILATIONPATH) 
UI_DIR = $(YADECOMPILATIONPATH) 
OBJECTS_DIR = $(YADECOMPILATIONPATH) 
QMAKE_LIBDIR = ../../../../../../plugins/Data/Body/InteractingGeometry/PolyhedralSweptSphere/$(YADEDYNLIBPATH) \
               ../../../../../../plugins/Data/Interaction/NarrowInteractionGeometry/MacroMicroContactGeometry/$(YADEDYNLIBPATH) \
               /usr/local/lib/yade/yade-libs/ \
               $(YADEDYNLIBPATH) 
QMAKE_CXXFLAGS_RELEASE += -lpthread \
                          -pthread 
QMAKE_CXXFLAGS_DEBUG += -lpthread \
                        -pthread 
DESTDIR = $(YADEDYNLIBPATH) 
CONFIG += debug \
          warn_on \
          dll 
TEMPLATE = lib 
HEADERS += SwiftPolyhedronProximityModeler.hpp 
SOURCES += SwiftPolyhedronProximityModeler.cpp 
