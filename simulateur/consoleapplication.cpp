#include "consoleapplication.h"

#include <spdlog/spdlog.h>

#include "magneticsimulation.h"

ConsoleApplication::ConsoleApplication(ConsoleApplicationParameters parameters, QObject* parent) :
	QObject(parent),
	m_parameters(std::move(parameters))
{
	spdlog::set_level(spdlog::level::trace);
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
	MagneticSimulation sim(m_parameters.inputFile);

	if (sim.checkValid())
	{
		sim.runSimulation();

		sim.saveResults("matriceA.png",
		                "matriceBr.png",
		                "matriceBz.png",
		                "matriceB.png",
		                "matriceAcontour.png");
	}

	return true;
}
