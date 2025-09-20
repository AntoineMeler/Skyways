QT             += gui
CONFIG         += c++17 console optimize_full
CONFIG         -= app_bundle
OBJECTS_DIR     = obj
QMAKE_CXXFLAGS += -fopenmp
DEFINES        += QT_DEPRECATED_WARNINGS
SOURCES        += main.cpp data.cpp particles.cpp image.cpp color.cpp mymath.cpp export_ply.cpp linking.cpp
LIBS           += -fopenmp -larmadillo

INCLUDEPATH += /opt/local/include /usr/include/eigen3

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
