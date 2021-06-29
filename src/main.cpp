#include <cstdlib>
#include <iostream>
#include <QApplication>
#include <QOpenGLContext>
#include "Display.h"
#include "commit.h"

int main(int argc, char * argv[])
{
	std::cout << "Qt version: " << qVersion() << std::endl;
	
	QSurfaceFormat fmt;
	if (QOpenGLContext::openGLModuleType() == QOpenGLContext::LibGL) 
	{
		std::cout << "OpenGL 3.3 context" << std::endl;
		fmt.setVersion(3, 3);
		fmt.setProfile(QSurfaceFormat::CoreProfile);
	}
	else 
	{
		std::cout << "OpenGL 3.0 context" << std::endl;
		fmt.setVersion(3, 0);
	}

	std::cout << "OpenGL Version: " << fmt.version().first << "." <<
	fmt.version().second << std::endl;
	QSurfaceFormat::setDefaultFormat(fmt);

	QCoreApplication::setAttribute(Qt::AA_ShareOpenGLContexts, true);
	QApplication app(argc, argv);

	QOpenGLContext *global = QOpenGLContext::globalShareContext();
	global->setShareContext(global);
	global->create();

	setlocale(LC_NUMERIC, "C");
	srand(time(NULL));
	
	std::cout << "Program version: " << CHECK_VERSION_COMMIT_ID << std::endl;
	
	Display display(argc, argv);
	display.show();
	/*
	SerumView serum(NULL);
	serum.show();
	serum.setCommandLineArgs(argc, argv);
	*/

	int status = app.exec();
	
	return status;
}
