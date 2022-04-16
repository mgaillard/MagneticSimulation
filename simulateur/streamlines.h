#pragma once

#include <Eigen/Core>

#include "scene.h"

/**
 * \brief Adapter for sampling the magnetic field as a vector field
 */
class MagneticVectorFieldAdapter
{
public:
    MagneticVectorFieldAdapter(Scene scene, Eigen::MatrixXd matBr, Eigen::MatrixXd matBz);

    /**
     * \brief Evaluate the vector field at position p and time t
     * \param t Time in seconds (constant in this case)
     * \param p 2D position in meters
     * \return The 2D vector sampled from the vector field
     */
    Eigen::Vector2d operator()(double t, const Eigen::Vector2d& p) const;

    /**
     * \brief Check if the point p is within bounds of the vector field
     * \param p 2D position in meters
     * \return True if the position is within bounds, or false
     */
    bool isWithinBounds(const Eigen::Vector2d& p) const;

private:
    const Scene m_scene;

    const Eigen::MatrixXd m_matBr;
    const Eigen::MatrixXd m_matBz;

    const int m_rows;
    const int m_cols;
};

/**
 * \brief Helper class to integrate a vector field using Runge-Kutta fourth order
 */
class IntegratorRK4
{
public:

    using VectorField = std::function<Eigen::Vector2d(double, const Eigen::Vector2d&)>;

    explicit IntegratorRK4(VectorField f, double h);

    Eigen::Vector2d evaluate(double t, const Eigen::Vector2d& p) const;

private:

    const VectorField m_func;
    const double m_h;
};

class StreamlinesPlacement
{
public:
    using Streamline = std::vector<Eigen::Vector2d>;

    explicit StreamlinesPlacement(Scene scene, Eigen::MatrixXd matBr, Eigen::MatrixXd matBz);

    int nbStreamlines() const;

    const Streamline& streamline(int index) const;

    void placeStreamline(const Eigen::Vector2d& seed);

private:
    MagneticVectorFieldAdapter m_vectorFieldAdapter;

    std::vector<Streamline> m_streamlines;
};
