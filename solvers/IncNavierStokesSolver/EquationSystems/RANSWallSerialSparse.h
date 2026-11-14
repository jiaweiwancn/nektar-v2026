///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_SERIAL_SPARSE_H
#define NEKTAR_RANS_WALL_SERIAL_SPARSE_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.h>
#include <LibUtilities/Communication/Comm.h>
#include <cstdint>
#include <memory>

namespace Nektar::RANSWallAlgebra
{
/** Deterministic row assembly -> immutable CSR. Never stages a global dense
 * matrix. Local contribution order is the dense assembler's original order.
 * Copies of finalized matrices share immutable storage, including mass caches.
 */
class SparseMatrix
{
public:
    SparseMatrix()=default;
    SparseMatrix(int rows,int columns):m_rows(rows),m_columns(columns)
    { Require(rows>=0 && columns>=0,"negative sparse shape");m_assembly.resize(rows); }
    int Rows() const { return m_rows; }
    int Columns() const { return m_columns; }
    double &Entry(int row,int column);
    void ClearRow(int row);
    double Get(int row,int column) const;
    void Finalize() const;
    std::size_t NonZeros() const { Finalize();return m_data->values.size(); }
    template<class Function> void Each(int row,Function function) const
    {
        Finalize();Require(row>=0 && row<m_rows,"sparse row index");
        for(int k=m_data->ptr[row];k<m_data->ptr[row+1];++k)
            function(m_data->ids[k],m_data->values[k]);
    }
    Vector Apply(const Vector &x) const;
    SparseMatrix Block(int row0,int row1,int col0,int col1) const;
    SparseMatrix Restrict(const std::vector<int> &ids) const;
    DenseMatrix DenseForVerification(int limit=1000) const;
    void Write(const std::string &file) const;
    const std::vector<int> &Pointers() const { Finalize();return m_data->ptr; }
    const std::vector<int> &Indices() const { Finalize();return m_data->ids; }
    const Vector &Values() const { Finalize();return m_data->values; }
private:
    struct Data { std::vector<int> ptr,ids;Vector values; };
    int m_rows=0,m_columns=0;
    mutable std::vector<std::map<int,double>> m_assembly;
    mutable std::shared_ptr<const Data> m_data;
};

class SparseFactor
{
public:
    explicit SparseFactor(const SparseMatrix &matrix,bool positiveMass=false);
    Vector Apply(const Vector &rhs) const;
    std::size_t FactorNonZeros() const;
private:
    struct Impl;
    std::shared_ptr<Impl> m_impl;
};

struct SparseKrylovResult { Vector solution;int iterations=0;double trueResidual=0; };
SparseKrylovResult SparsePressureGMRES(const LibUtilities::SessionReaderSharedPtr &session,
    const LibUtilities::CommSharedPtr &comm,const Vector &rhs,
    const std::function<Vector(const Vector&)> &action,
    const std::function<Vector(const Vector&)> &preconditioner);
}
#endif
