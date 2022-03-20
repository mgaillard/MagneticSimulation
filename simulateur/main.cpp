#include <QCoreApplication>
#include <QCommandLineParser>
#include <QFileInfo>
#include <QDebug>

#include "consoleapplication.h"

enum class CommandLineParseResult
{
    Ok,
    Error
};

CommandLineParseResult parseCommandLine(const QCoreApplication& app,
                                        QCommandLineParser& parser,
                                        ConsoleApplicationParameters* parameters);

#include <Eigen/Dense>
#include <iostream>
int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);
    
    // Parse command line arguments
    QCommandLineParser parser;
    parser.setApplicationDescription("Magnetic Simulation");
    ConsoleApplicationParameters parameters;
    const CommandLineParseResult action = parseCommandLine(app, parser, &parameters);

    if (action == CommandLineParseResult::Ok)
    {
        ConsoleApplication consoleApp(parameters);

        // Use a queued connection because the console app
        // may emit the finished signal before app.exec() is called
        QObject::connect(&consoleApp, &ConsoleApplication::finished,
                         &app, &QCoreApplication::exit, Qt::QueuedConnection);

        consoleApp.exec();
        return app.exec();
    }
    
    return EXIT_FAILURE;
}

CommandLineParseResult parseCommandLine(
	const QCoreApplication& app,
	QCommandLineParser& parser,
	ConsoleApplicationParameters* parameters)
{
	const QCommandLineOption helpOption = parser.addHelpOption();


	// An option to set the input file
	const QCommandLineOption inputOption(
		QStringList() << "i" << "input",
		QCoreApplication::translate("main", "Input file."),
		QCoreApplication::translate("main", "file"));
	parser.addOption(inputOption);

	// Process the actual command line arguments given by the user
	parser.process(app);

	// No input file given, this is an error
	if (!parser.isSet(inputOption))
	{
		qWarning() << "Error: no input file was given.";
		return CommandLineParseResult::Error;
	}

	// Input file given
	if (parser.isSet(inputOption))
	{
		// Check if the input file exists
		const QFileInfo fileInfo(parser.value(inputOption));

		if (!fileInfo.exists())
		{
			qWarning() << "Error: the input file does not exist.";
			return CommandLineParseResult::Error;
		}
		parameters->inputFile = fileInfo.absoluteFilePath();
		return CommandLineParseResult::Ok;
	}

	return CommandLineParseResult::Error;
}
