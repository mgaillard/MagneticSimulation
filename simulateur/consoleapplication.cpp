#include "consoleapplication.h"

#include <spdlog/spdlog.h>

#include "simulateur.h"

ConsoleApplication::ConsoleApplication(ConsoleApplicationParameters parameters, QObject* parent) :
	QObject(parent),
	m_parameters(std::move(parameters))
{

}

void ConsoleApplication::exec()
{
	// returnCode is either 0 if the output is successfully saved, or 1 if a problem occurred
	int returnCode = 1;

	if (run())
	{
		returnCode = 0;
	}

	emit finished(returnCode);
}

bool ConsoleApplication::run() const
{
	Simulateur sim(m_parameters.inputFile);

	if (sim.validationSimulation())
	{
		sim.simuler();

		sim.enregistrerResultats(QString("matriceA.png"),
		                         QString("matriceBr.png"),
		                         QString("matriceBz.png"),
		                         QString("matriceB.png"),
								 QString("matriceAcontour.png"));
	}

	return true;
}
