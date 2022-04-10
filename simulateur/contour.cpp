#include "contour.h"

#include <bitset>

// ---------------------------------------- Private functions ----------------------------------------

double lerp(double a, double b, double x)
{
	return a + x * (b - a);
}

Eigen::Vector2d lerp(const Eigen::Vector2i& a, const Eigen::Vector2i& b, double x)
{
	return {
		lerp(a.x(), b.x(), x),
		lerp(a.y(), b.y(), x)
	};
}

Eigen::Vector2d lerpIsoValue(
	const Eigen::MatrixXd& matrix,
	const Eigen::Vector2i& a,
	const Eigen::Vector2i& b,
	double k)
{
	const double az = matrix(a.x(), a.y());
	const double bz = matrix(b.x(), b.y());

	if ((az <= k && k <= bz) || (az >= k && k >= bz))
	{
		const auto x = (k - az) / (bz - az);

		return lerp(a, b, x);
	}

	return lerp(a, b, 0.5);
}

// ---------------------------------------- Public functions ----------------------------------------

void Contour::addContour(const Eigen::MatrixXd& matrix, double k)
{
	for (int i = 0; i < matrix.rows() - 1; i++)
	{
		for (int j = 0; j < matrix.cols() - 1; j++)
		{
			const Eigen::Vector2i blPoint(i + 1, j);
			const Eigen::Vector2i brPoint(i + 1, j + 1);
			const Eigen::Vector2i tlPoint(i, j);
			const Eigen::Vector2i trPoint(i, j + 1);

			// Pack the configuration in 4 bits
			std::bitset<4> config;
			config.set(0, matrix(blPoint.x(), blPoint.y()) >= k);
			config.set(1, matrix(brPoint.x(), brPoint.y()) >= k);
			config.set(2, matrix(tlPoint.x(), tlPoint.y()) >= k);
			config.set(3, matrix(trPoint.x(), trPoint.y()) >= k);

			// Mirror configurations 
			if (config.to_ulong() > 7)
			{
				config.flip();
			}

			switch (config.to_ulong())
			{
			case 0:
				// All corners have the same value: no contour
				break;
			case 1:
				{
					// Bottom left activated, the rest deactivated
					// Edge between left and bottom
					const auto left = lerpIsoValue(matrix, blPoint, tlPoint, k);
					const auto bottom = lerpIsoValue(matrix, blPoint, brPoint, k);
					addLine(left, bottom);
				}
				break;
			case 2:
				{
					// Bottom right activated, the rest deactivated
					// Edge between bottom and right
					const auto bottom = lerpIsoValue(matrix, blPoint, brPoint, k);
					const auto right = lerpIsoValue(matrix, brPoint, trPoint, k);
					addLine(bottom, right);
				}
				break;
			case 3:
				{
					// Bottom activated, top deactivated
					// Edge between left and right
					const auto left = lerpIsoValue(matrix, blPoint, tlPoint, k);
					const auto right = lerpIsoValue(matrix, brPoint, trPoint, k);
					addLine(left, right);
				}
				break;
			case 4:
				{
					// Top left activated, the rest deactivated
					// Edge between left and top
					const auto left = lerpIsoValue(matrix, blPoint, tlPoint, k);
					const auto top = lerpIsoValue(matrix, tlPoint, trPoint, k);
					addLine(left, top);
				}
				break;
			case 5:
				{
					// Left activated, right deactivated
					// Edge between bottom and top
					const auto bottom = lerpIsoValue(matrix, blPoint, brPoint, k);
					const auto top = lerpIsoValue(matrix, tlPoint, trPoint, k);
					addLine(bottom, top);
				}
				break;
			case 6:
				{
					const double blPointZ = matrix(blPoint.x(), blPoint.y());
					const double brPointZ = matrix(brPoint.x(), brPoint.y());
					const double tlPointZ = matrix(tlPoint.x(), tlPoint.y());
					const double trPointZ = matrix(trPoint.x(), trPoint.y());

					// Ambiguous case: two corners activated, two corners deactivated
					// Value of z at the saddle point s
					const auto sz = (blPointZ * trPointZ + brPointZ * tlPointZ)
					                / (blPointZ - brPointZ + trPointZ - tlPointZ);

					const auto bottom = lerpIsoValue(matrix, blPoint, brPoint, k);
					const auto right = lerpIsoValue(matrix, brPoint, trPoint, k);
					const auto left = lerpIsoValue(matrix, blPoint, tlPoint, k);
					const auto top = lerpIsoValue(matrix, tlPoint, trPoint, k);

					// The saddle point has the same sign as the bottom left corner
					if (sz >= k && blPointZ >= k)
					{
						addLine(bottom, right);
						addLine(left, top);
					}
					// The saddle point has the same sign as the bottom right corner
					else if (sz < k && blPointZ < k)
					{
						addLine(bottom, left);
						addLine(right, top);
					}
					else
					{
						// No luck, the case is really ambiguous
					}
				}
				break;
			case 7:
				{
					// Top right activated, the rest deactivated
					// Edge between top and right
					const auto top = lerpIsoValue(matrix, tlPoint, trPoint, k);
					const auto right = lerpIsoValue(matrix, brPoint, trPoint, k);
					addLine(top, right);
				}
				break;
			default:
				break;
			}
		}
	}
}

const std::vector<Eigen::Vector2d>& Contour::points() const
{
	return m_points;
}

const std::vector<std::vector<int>>& Contour::curves() const
{
	return m_curves;
}

void Contour::addLine(const Eigen::Vector2d& a, const Eigen::Vector2d& b)
{
	const auto [newA, indexA] = addPoint(a);
	const auto [newB, indexB] = addPoint(b);

	// If a is new, b is new, new curve
	if (newA && newB)
	{
		m_curves.push_back({ indexA, indexB });
	}
	// If a is not new, b is new, find the curve with a in it
	else if (!newA && newB)
	{
		const int curveIndex = findCurveWithPoint(indexA);
		assert(curveIndex >= 0);

		// a is either at the beginning or the end of the curve
		if (atBeginningOfCurve(indexA, curveIndex))
		{
			// Insert b
			m_curves[curveIndex].insert(m_curves[curveIndex].begin(), 1, indexB);
		}
		else if (atEndOfCurve(indexA, curveIndex))
		{
			// Append b
			m_curves[curveIndex].push_back(indexB);
		}
		else
		{
			// The point a should not be in the middle of the curve, that means there is a T junction
			assert(false);
		}
	}
	// If a is new, b is not new, find the curve with b in it
	else if (newA && !newB)
	{
		const int curveIndex = findCurveWithPoint(indexB);
		assert(curveIndex >= 0);

		// b is either at the beginning or the end of the curve
		if (atBeginningOfCurve(indexB, curveIndex))
		{
			// Insert a
			m_curves[curveIndex].insert(m_curves[curveIndex].begin(), 1, indexA);
		}
		else if (atEndOfCurve(indexB, curveIndex))
		{
			// Append a
			m_curves[curveIndex].push_back(indexA);
		}
		else
		{
			// The point a should not be in the middle of the curve, that means there is a T junction
			assert(false);
		}
	}
	// If both a and b are not new, find the curves with a and b in them, and merge
	else
	{
		const int curveIndexA = findCurveWithPoint(indexA);
		const int curveIndexB = findCurveWithPoint(indexB);
		assert(curveIndexA >= 0 && curveIndexB >= 0);

		// If both points are on the same curve
		if (curveIndexA == curveIndexB)
		{
			// a at the end and b at the beginning
			if (atEndOfCurve(indexA, curveIndexA)
				&& atBeginningOfCurve(indexB, curveIndexB))
			{
				// Append b
				m_curves[curveIndexB].push_back(indexB);
			}
			// a at the beginning and b at the end
			else if (atBeginningOfCurve(indexA, curveIndexA)
				&& atEndOfCurve(indexB, curveIndexB))
			{
				// Append a
				m_curves[curveIndexA].push_back(indexA);
			}
			else
			{
				// The points should not be in the middle of a curve, that means there is a T junction
				assert(false);
			}
		}
		else
		{
			// Both points at the beginning of the curve
			if (atBeginningOfCurve(indexA, curveIndexA)
				&& atEndOfCurve(indexB, curveIndexB))
			{
				// Append the reverse of curve a to curve b
				m_curves[curveIndexB].insert(m_curves[curveIndexB].end(),
					m_curves[curveIndexA].rbegin(),
					m_curves[curveIndexA].rend());

				// Delete curve a
				m_curves.erase(m_curves.begin() + curveIndexA);
			}
			// a at the end and b at the beginning
			else if (atEndOfCurve(indexA, curveIndexA)
				&& atBeginningOfCurve(indexB, curveIndexB))
			{
				// Append curve b at the end of curve a
				m_curves[curveIndexA].insert(m_curves[curveIndexA].end(),
					m_curves[curveIndexB].begin(),
					m_curves[curveIndexB].end());

				// Delete curve b
				m_curves.erase(m_curves.begin() + curveIndexB);
			}
			// a at the beginning and b at the end
			else if (atBeginningOfCurve(indexA, curveIndexA)
				&& atEndOfCurve(indexB, curveIndexB))
			{
				// Append curve a at the end of curve b
				m_curves[curveIndexB].insert(m_curves[curveIndexB].end(),
					m_curves[curveIndexA].begin(),
					m_curves[curveIndexA].end());

				// Delete curve a
				m_curves.erase(m_curves.begin() + curveIndexA);
			}
			// both points at the end of the curve
			else if (atEndOfCurve(indexA, curveIndexA)
				&& atEndOfCurve(indexB, curveIndexB))
			{
				// Append the reverse of curve b to curve a
				m_curves[curveIndexA].insert(m_curves[curveIndexA].end(),
					m_curves[curveIndexB].rbegin(),
					m_curves[curveIndexB].rend());

				// Delete curve b
				m_curves.erase(m_curves.begin() + curveIndexB);
			}
			else
			{
				// The points should not be in the middle of a curve, that means there is a T junction
				assert(false);
			}
		}
	}
}

std::tuple<bool, int> Contour::addPoint(const Eigen::Vector2d& p)
{
	constexpr auto eps = 1e-12;
	
	for (unsigned int i = 0; i < m_points.size(); i++)
	{
		if ((m_points[i] - p).norm() <= eps)
		{
			return {false, static_cast<int>(i)};
		}
	}

	const int index = static_cast<int>(m_points.size());
	m_points.push_back(p);
	return {true, index};
}

int Contour::findCurveWithPoint(int index) const
{
	for (unsigned int i = 0; i < m_curves.size(); i++)
	{
		if (std::find(m_curves[i].begin(), m_curves[i].end(), index) != m_curves[i].end())
		{
			return static_cast<int>(i);
		}
	}

	return -1;
}

bool Contour::atBeginningOfCurve(int pointIndex, int curveIndex) const
{
	return m_curves[curveIndex].front() == pointIndex;
}

bool Contour::atEndOfCurve(int pointIndex, int curveIndex) const
{
	return m_curves[curveIndex].back() == pointIndex;
}
