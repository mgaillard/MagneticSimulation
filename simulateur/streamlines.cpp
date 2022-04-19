#include "streamlines.h"

#include "spdlog/spdlog.h"

MagneticVectorFieldAdapter::MagneticVectorFieldAdapter(
	Scene scene,
	Eigen::MatrixXd matBr,
	Eigen::MatrixXd matBz) :
	m_scene(std::move(scene)),
	m_matBr(std::move(matBr)),
	m_matBz(std::move(matBz)),
	m_rows(static_cast<int>(std::min(m_matBr.rows(), m_matBz.rows()))),
	m_cols(static_cast<int>(std::min(m_matBr.cols(), m_matBz.cols())))
{

}

Eigen::Vector2d MagneticVectorFieldAdapter::operator()(
	[[maybe_unused]] double t,
	const Eigen::Vector2d& p) const
{
	const auto gridCoords = m_scene.convertToGridCoords(p);

	const auto x = static_cast<int>(gridCoords.x());
	const auto y = static_cast<int>(gridCoords.y());

	if (x >= 0 && x < m_cols && y >= 0 && y < m_rows)
	{
		// TODO: use bi-cubic sampling
		// Sample nearest point
		return {
			m_matBr(y, x),
			m_matBz(y, x),
		};
	}

	return { 0.0, 0.0 };
}

bool MagneticVectorFieldAdapter::isWithinBounds(const Eigen::Vector2d& p) const
{
	return m_scene.isInside(p);
}

IntegratorRK4::IntegratorRK4(VectorField f, double h): m_func(std::move(f)), m_h(h)
{
	
}

Eigen::Vector2d IntegratorRK4::evaluate(double t, const Eigen::Vector2d& p) const
{
	const Eigen::Vector2d k1 = m_func(t, p);
	const Eigen::Vector2d k2 = m_func(t + m_h / 2.0, p + (m_h / 2.0) * k1);
	const Eigen::Vector2d k3 = m_func(t + m_h / 2.0, p + (m_h / 2.0) * k2);
	const Eigen::Vector2d k4 = m_func(t + m_h, p + m_h * k3);

	return p + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * (m_h / 6.0);
}

StreamlinesPlacement::StreamlinesPlacement(Scene scene, Eigen::MatrixXd matBr, Eigen::MatrixXd matBz) :
	m_vectorFieldAdapter(std::move(scene), std::move(matBr), std::move(matBz))
{

}

int StreamlinesPlacement::nbStreamlines() const
{
	return static_cast<int>(m_streamlines.size());
}

const StreamlinesPlacement::Streamline& StreamlinesPlacement::streamline(int index) const
{
	return m_streamlines.at(index);
}

void StreamlinesPlacement::placeStreamline(const Eigen::Vector2d& seed)
{
	// Maximum number of integration steps
	constexpr int maxNumberSteps = 100000;
	// Minimum distance to the current streamline below which integration is cancelled
	constexpr double minDistance = 1e-3;

	double t = 0.0;
	constexpr double h = 1.0;
	const IntegratorRK4 backwardIntegrator(m_vectorFieldAdapter, -h);
	const IntegratorRK4 forwardIntegrator(m_vectorFieldAdapter, h);

	Streamline streamline(1, seed);
	
	// Backward integration
	auto point = streamline.back();
	for (int i = 0; i < maxNumberSteps; i++)
	{
		point = backwardIntegrator.evaluate(t, point);
		t -= h;

		// If we go outside the scene, stop
		if (!m_vectorFieldAdapter.isWithinBounds(point))
		{
			break;
		}
		
		// Distance to the beginning of the streamline
		const auto distToBeginning = (streamline.front() - point).norm();

		// If too close to the other end of the streamline, stop
		if (distToBeginning < minDistance / 2.0)
		{
			// TODO: early stop the streamline when it loops on itself
			// break;
		}

		// Distance to the previous point in the streamline
		const auto distToPrevious = (streamline.back() - point).norm();

		// Only add if there is a certain distance from the last point
		// This prevent the streamline to have a resolution that is too fine
		if (distToPrevious >= minDistance)
		{
			streamline.emplace_back(point);
		}
	}

	// Reverse streamline
	std::reverse(streamline.begin(), streamline.end());

	// Forward integration
	point = streamline.back();
	for (int i = 0; i < maxNumberSteps; i++)
	{
		point = forwardIntegrator.evaluate(t, point);
		t += h;

		// If we go outside the scene, stop
		if (!m_vectorFieldAdapter.isWithinBounds(point))
		{
			break;
		}

		// Distance to the beginning of the streamline
		const auto distToBeginning = (streamline.front() - point).norm();

		// If too close to the other end of the streamline, stop
		if (distToBeginning < minDistance / 2.0)
		{
			// TODO: early stop the streamline when it loops on itself
			// break;
		}

		// Distance to the previous point in the streamline
		const auto distToPrevious = (streamline.back() - point).norm();

		// Only add if there is a certain distance from the last point
		// This prevent the streamline to have a resolution that is too fine
		if (distToPrevious >= minDistance)
		{
			streamline.emplace_back(point);
		}
	}

	m_streamlines.push_back(streamline);
}
