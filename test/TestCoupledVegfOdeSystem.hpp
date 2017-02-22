/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTCOUPLEDVEGFODESYSTEM_HPP_
#define TESTCOUPLEDVEGFODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractFeVolumeIntegralAssembler.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "MassMatrixAssembler.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"
#include "ConstBoundaryCondition.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "RobinConditionsSurfaceTermAssembler.hpp"
#include "AbstractFeSurfaceIntegralAssemblerWithMatrix.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "Debug.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * A simple parabolic PDE used in tests.
 */
template <int SPACE_DIM>
class HeatEquation : public AbstractLinearParabolicPde<SPACE_DIM>
{

public:
    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& , double, Element<SPACE_DIM,SPACE_DIM>* pElement=NULL)
    {
        return 0.0;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& , Element<SPACE_DIM,SPACE_DIM>* pElement=NULL)
    {
        return identity_matrix<double>(SPACE_DIM);
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& )
    {
        return 1;
    }
};

template<unsigned DIM>
class SurfaceIntegralCalculator
    : public AbstractFeSurfaceIntegralAssemblerWithMatrix<DIM,DIM,1>
{

public:

    double CalculateSurfaceIntegral()
    {
        double integral_total = 0.0;

        typename BoundaryConditionsContainer<DIM,DIM,1>::NeumannMapIterator
            neumann_iterator = this->mpBoundaryConditions->BeginNeumann();

        // Iterate over defined conditions
        while (neumann_iterator != this->mpBoundaryConditions->EndNeumann())
        {
            const BoundaryElement<DIM-1,DIM>& r_surf_element = *(neumann_iterator->first);

            c_vector<double, DIM> weighted_direction;
            double jacobian_determinant;
            this->mpMesh->GetWeightedDirectionForBoundaryElement(r_surf_element.GetIndex(), weighted_direction,
                    jacobian_determinant);

            // Allocate memory for the basis function values
            c_vector<double, DIM> phi;

            // Loop over Gauss points
            double element_integral = 0.0;
            for (unsigned quad_index=0; quad_index<this->mpSurfaceQuadRule->GetNumQuadPoints(); quad_index++)
            {
                const ChastePoint<DIM-1>& quad_point = this->mpSurfaceQuadRule->rGetQuadPoint(quad_index);

                LinearBasisFunction<DIM-1>::ComputeBasisFunctions(quad_point, phi);
                double node_contribution = 0.0;
                for (unsigned i=0; i<r_surf_element.GetNumNodes(); i++)
                {
                    double u_at_node = this->GetCurrentSolutionOrGuessValue(r_surf_element.GetNode(i)->GetIndex(), 0);
                    node_contribution += phi(i) * u_at_node;
                }
                double wJ = jacobian_determinant * this->mpSurfaceQuadRule->GetWeight(quad_index);
                element_integral += wJ*node_contribution;
            }
            integral_total += element_integral;
            ++neumann_iterator;
        }
        return integral_total;
    }

public:

    SurfaceIntegralCalculator(AbstractTetrahedralMesh<DIM,DIM>* pMesh,
            BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions)
        : AbstractFeSurfaceIntegralAssemblerWithMatrix<DIM,DIM,1>(pMesh, pBoundaryConditions)
    {

    }
};

template<unsigned DIM>
class ParabolicTermsAssembler
    : public AbstractFeVolumeIntegralAssembler<DIM,DIM,1 ,true ,true, NORMAL>
{

protected:

    AbstractLinearParabolicPde<DIM,DIM>* mpParabolicPde;

private:

    /* Even when a class isn't being written for a very general dimensions sometimes it is a good idea
     * to define the following, and then use `ELEMENT_DIM` etc in the below, as it can make the code a
     * bit easier to understand.
     */
    static const unsigned ELEMENT_DIM = DIM;
    static const unsigned SPACE_DIM = DIM;
    static const unsigned PROBLEM_DIM = 1;

    /* We are assembling a matrix, we means we need to provide a `ComputeMatrixTerm()` method, to return the
     * elemental contribution to the RHS matrix. Note that `ELEMENT_DIM+1` is the number of
     * nodes in the element (=number of basis functions).
     */
    c_matrix<double,PROBLEM_DIM*(ELEMENT_DIM+1),PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(
                                                                                c_vector<double, ELEMENT_DIM+1> &rPhi,
                                                                                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                                                                                ChastePoint<SPACE_DIM> &rX,
                                                                                c_vector<double,PROBLEM_DIM> &rU,
                                                                                c_matrix<double, PROBLEM_DIM, SPACE_DIM> &rGradU /* not used */,
                                                                                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        c_matrix<double, SPACE_DIM, SPACE_DIM> pde_diffusion_term = mpParabolicPde->ComputeDiffusionTerm(rX, pElement);

        return    prod( trans(rGradPhi), c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
                + PdeSimulationTime::GetPdeTimeStepInverse() * mpParabolicPde->ComputeDuDtCoefficientFunction(rX) * outer_prod(rPhi, rPhi);

    }

    /**
     * @return the term to be added to the element stiffness vector - see AbstractFeVolumeIntegralAssembler
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)

        {
            return (mpParabolicPde->ComputeSourceTerm(rX, rU(0), pElement)
                    + PdeSimulationTime::GetPdeTimeStepInverse() * mpParabolicPde->ComputeDuDtCoefficientFunction(rX) * rU(0)) * rPhi;
        }

public:
    ParabolicTermsAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh,
            AbstractLinearParabolicPde<DIM,DIM>* pPde)
        : AbstractFeVolumeIntegralAssembler<DIM,DIM,1,true,true,NORMAL>(pMesh),
          mpParabolicPde(pPde)
    {
    }
};

template<unsigned DIM>
class CoupledVegfOdeSystemSolver : public AbstractDynamicLinearPdeSolver<DIM,DIM,1>
{

protected:

    AbstractLinearParabolicPde<DIM,DIM>* mpParabolicPde;

private:

    /* The constuctor will take in a mesh and a BCC, the latter will be stored as a member variable */
    BoundaryConditionsContainer<DIM,DIM,1>* mpBoundaryConditions;

    Mat mRhsRobinMatrix;

    double mPermeability;

    /* This is the main method which needs to be implemented. It takes in the current solution, and a
     * boolean saying whether the matrix (ie A in Ax=b) is being computed or not.
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {

        // Update the surface integral value
        SurfaceIntegralCalculator<DIM> surf_calc(this->mpMesh, mpBoundaryConditions);
        surf_calc.SetCurrentSolution(currentSolution);
        double integral_value = surf_calc.CalculateSurfaceIntegral();
        std::cout << integral_value << std::endl;

        // Update the Boundary condition value



        /* This is how to use assemblers to set up matrices. We declare a mass matrix assembler,
         * pass it the LHS matrix of the linear system, and tell it to assemble. We also declare
         * one of our purpose-built `RhsMatrixAssemblers`, pass it the matrix `mRhsMatrix`, and
         * tell it to assemble.
         *
         * '''Important note''': if any of the assemblers will require the current solution (ie solution
         * at the current timestep), this needs to be passed to the assembler, as in the commented
         * line below.
         */
        ParabolicTermsAssembler<DIM> parabolic_terms_assembler(this->mpMesh, this->mpParabolicPde);

        if (computeMatrix)
        {
            parabolic_terms_assembler.SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
            parabolic_terms_assembler.AssembleMatrix();
            this->mpLinearSystem->FinaliseLhsMatrix(); // (Petsc communication)
        }
        else
        {
            RobinConditionsSurfaceTermAssembler<DIM,DIM,1> surface_integral_assembler(this->mpMesh, mpBoundaryConditions);
            surface_integral_assembler.SetMatrixToAssemble(mRhsRobinMatrix, true);
            surface_integral_assembler.AssembleMatrix();
            PetscMatTools::Finalise(mRhsRobinMatrix);
            MatMult(mRhsRobinMatrix, currentSolution, this->mpLinearSystem->rGetRhsVector());

            surface_integral_assembler.SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false);
            surface_integral_assembler.AssembleVector();
            this->mpLinearSystem->FinaliseRhsVector(); // (Petsc communication)

            parabolic_terms_assembler.SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false);
            parabolic_terms_assembler.SetCurrentSolution(currentSolution);
            parabolic_terms_assembler.AssembleVector();
            this->mpLinearSystem->FinaliseRhsVector(); // (Petsc communication)
        }


//        /* The third assembler we use is the `NaturalNeumannSurfaceTermAssembler`, which assembles
//         * the vector `c` defined above, using the Neumann BCs stored in the `BoundaryConditionsContainer`
//         * which is passed in in the constructor
//         */
//        NaturalNeumannSurfaceTermAssembler<DIM,DIM,1> surface_integral_assembler(this->mpMesh, mpBoundaryConditions);
//        surface_integral_assembler.SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false /*don't zero vector before assembling!*/);
//        surface_integral_assembler.Assemble();

        /* Some necessary PETSc communication before applying Dirichet BCs */
        this->mpLinearSystem->FinaliseRhsVector();         // (Petsc communication)
        this->mpLinearSystem->SwitchWriteModeLhsMatrix();  // (Petsc communication - needs to called when going from adding entries to inserting entries)

        /* Apply the dirichlet BCs from the BCC to the linear system */
        mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

        /* Some necessary PETSc communication to finish */
        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->FinaliseLhsMatrix();

        // Update the Pellet VEGF for the next step
    }

public:
    /* The constructor needs to call the parent constructor, save the BCC, ''say that the (LHS) matrix is constant
     * in time'' (so it is only computed once), and allocate memory for the RHS matrix.
     */
    CoupledVegfOdeSystemSolver(TetrahedralMesh<DIM,DIM>* pMesh,
                               BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions,
                               AbstractLinearParabolicPde<DIM,DIM>* pPde)
         : AbstractDynamicLinearPdeSolver<DIM,DIM,1>(pMesh),
           mpParabolicPde(pPde),
           mpBoundaryConditions(pBoundaryConditions),
           mPermeability(1.0)
    {
        this->mMatrixIsConstant = true;
        PetscTools::SetupMat(mRhsRobinMatrix, this->mpMesh->GetNumNodes(), this->mpMesh->GetNumNodes(), 9);
    }

    ~CoupledVegfOdeSystemSolver()
    {
        PetscTools::Destroy(mRhsRobinMatrix);
    }

};

class TestCoupledVegfOdeSystem : public CxxTest::TestSuite
{
public:
    void TestExplicitSolver() throw (Exception)
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(2.0 /*h*/, 50.0 /*width*/, 18.0 /*height*/);

        HeatEquation<2> pde;

        // Set up BCs u=0 on entire boundary
        BoundaryConditionsContainer<2,2,1> bcc;

        // Hijack the Neumann boundary condition flags
        TetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* p_left_wall_boundary_condition = new ConstBoundaryCondition<2>(1.0);
        while (surf_iter != mesh.GetBoundaryElementIteratorEnd())
        {
            unsigned node_index = (*surf_iter)->GetNodeGlobalIndex(0);
            double x = mesh.GetNode(node_index)->GetPoint()[0];
            double y = mesh.GetNode(node_index)->GetPoint()[1];
            if (x < 0.2 and y>0.0)
            {
                bcc.AddNeumannBoundaryCondition(*surf_iter, p_left_wall_boundary_condition);
            }
            surf_iter++;
        }

        CoupledVegfOdeSystemSolver<2> solver(&mesh, &bcc, &pde);

        /* The interface is exactly the same as the `SimpleLinearParabolicSolver`. */
        solver.SetTimeStep(1.0);
        solver.SetTimes(0.0, 100.0);

        std::vector<double> init_cond(mesh.GetNumNodes(), 0.0);
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        solver.SetOutputDirectoryAndPrefix("CoupledVegfOdeSystem","results");
        solver.SetOutputToVtk(true);
        solver.SetPrintingTimestepMultiple(10);

        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Tidy up
        PetscTools::Destroy(initial_condition);
        PetscTools::Destroy(result);
    }
};

#endif // TESTCOUPLEDVEGFODESYSTEM_HPP_
