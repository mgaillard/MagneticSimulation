#include "solver.h"

Solver::Solver() : m_rows(0), m_cols(0)
{

}

void Solver::setSize(const int rows, const int cols)
{
    m_rows = rows;
    m_cols = cols;
}
