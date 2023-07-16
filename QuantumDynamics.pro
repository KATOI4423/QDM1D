QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17
CONFIG += c++14
CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
	../../qcustomplot/qcustomplot.cpp \
	main.cpp \
	mainwindow.cpp \
	quantum_dynamics_1d.cpp

HEADERS += \
	../../qcustomplot/qcustomplot.h \
	mainwindow.h \
	quantum_dynamics_1d.h

FORMS += \
	mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

QT       += printsupport
INCLUDEPATH += C:/Users/mrchi/Documents/MadeByQt/qcustomplot
INCLUDEPATH += C:/Users/mrchi/Documents/C++Header/eigen-3.4.0/Eigen

CONFIG += debug
CONFIG += debug_and_release
CONFIG += release

DESTDIR = build
win32 {
	TARGET_NAME = $${TARGET}.exe
}
else: macx {
	TARGET_NAME = $${TARGET}.app
}
else: linux {
	TARGET_NAME = $${TARGET}
}
TARGET_PATH = $$DESTDIR/$${TARGET_NAME}

win32 {
	QMAKE_POST_LINK =  $$dirname(QMAKE_QMAKE)/windeployqt $${TARGET_PATH} -qmldir=$$PWD --no-translations$$escape_expand(\n\t)
	QMAKE_POST_LINK += cd$$escape_expand(\n\t)
}
else: macx {
	QMAKE_POST_LINK =  $$dirname(QMAKE_QMAKE)/macdeployqt $${TARGET_PATH} -qmldir=$$PWD$$escape_expand(\n\t)
	QMAKE_POST_LINK += pwd$$escape_expand(\n\t)
}
else: linux {
	QMAKE_POST_LINK =  ~/Qt/linuxdeployqt-7-x86_64.AppImage $${TARGET_PATH} -qmldir=$$PWD$$escape_expand -unsupported-allow-new-glibc$$escape_expand(\n\t)
	QMAKE_POST_LINK += pwd$$escape_expand(\n\t)
}
