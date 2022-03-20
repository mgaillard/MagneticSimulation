#pragma once

#include <QObject>

struct ConsoleApplicationParameters
{
	QString inputFile;
};

class ConsoleApplication final : public QObject
{
	Q_OBJECT

public:
	explicit ConsoleApplication(ConsoleApplicationParameters parameters, QObject* parent = Q_NULLPTR);

signals:
	/**
	 * \brief This signal is emitted when the application has finished
	 * \param returnCode Return code of the application
	 */
	void finished(int returnCode);

public slots:
	/**
	 * \brief Execute the console application
	 */
	void exec();

private:
	/**
	 * \brief Run the application
	 * \return True if application completed successfully
	 */
	bool run() const;

	ConsoleApplicationParameters m_parameters;
};
