/******************************************************************************
*                    Anisotropy plugin for SOFA                               *
*                         version 1.0                                         *
*                       Copyright © Inria                                     *
*                       All rights reserved                                   *
*                       2024                                                  *
*                                                                             *
* This software is under the GNU General Public License v2 (GPLv2)            *
*            https://www.gnu.org/licenses/licenses.en.html                    *
*                                                                             *
*                                                                             *
*                                                                             *
* Authors: Felix Vanneste                                                     *
*                                                                             *
* Contact information: felix.vanneste@inria.information                       *
******************************************************************************/
#pragma once

#include <Anisotropy/component/initAnisotropy.h>

#include <sofa/core/topology/TopologyData.h>
//#include <sofa/core/topology/TopologyHandler.h>

#include <sofa/core/behavior/ForceField.h>

#include <sofa/type/MatSym.h>

namespace sofa::core::behavior
{

template< class T > class MechanicalState;

} // namespace sofa::core::behavior

namespace anisotropy::forcefield
{

using namespace sofa;
using namespace sofa::defaulttype;
using namespace sofa::core::topology;
using namespace sofa::type;

template<class DataTypes>
class TetrahedronAnisotropicForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TetrahedronAnisotropicForceField,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::Real        Real        ;
    typedef typename DataTypes::Coord       Coord       ;
    typedef typename DataTypes::Deriv       Deriv       ;
    typedef typename DataTypes::VecCoord    VecCoord    ;
    typedef typename DataTypes::VecDeriv    VecDeriv    ;
    typedef typename DataTypes::VecReal     VecReal     ;
    typedef Data<VecCoord>                  DataVecCoord;
    typedef Data<VecDeriv>                  DataVecDeriv;

    typedef core::topology::BaseMeshTopology::Tetrahedron Tetrahedron;
    typedef core::topology::BaseMeshTopology::TetraID TetraID;
    typedef core::topology::BaseMeshTopology::Tetra Tetra;
    typedef core::topology::BaseMeshTopology::Point Point;
    typedef core::topology::BaseMeshTopology::Triangle Triangle;
    typedef core::topology::BaseMeshTopology::Edge Edge;
    typedef core::topology::BaseMeshTopology::Quad Quad;
    typedef core::topology::BaseMeshTopology::Hexahedron Hexahedron;
    typedef core::topology::BaseMeshTopology::EdgesInTriangle EdgesInTriangle;
    typedef core::topology::BaseMeshTopology::EdgesInTetrahedron EdgesInTetrahedron;
    typedef core::topology::BaseMeshTopology::EdgesInQuad EdgesInQuad;
    typedef core::topology::BaseMeshTopology::EdgesInHexahedron EdgesInHexahedron;
    typedef core::topology::BaseMeshTopology::TrianglesInTetrahedron TrianglesInTetrahedron;
    typedef core::topology::BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;

    typedef type::vector<Coord> SetAnisotropyDirectionArray; // When the model is anisotropic, for instance in invariant I4
    typedef type::vector<Real> SetParameterArray; //necessary to store hyperelastic parameters (mu, lambda ...)

    typedef Mat<3,3,Real>       Mat3x3  ;
    typedef Mat<6,9,Real>       Mat6x9  ;


    typedef type::MatSym<3,Real> MatrixSym;
    // In case of non 3D template
    typedef Vec<3,Real> Vec3;
    typedef Vec<4,Real> Vec4;
    typedef Vec<5,Real> Vec5;
    typedef Vec<6,Real> Vec6;
    typedef Vec<9,Real> Vec9;
    typedef Vec<6,Mat3x3> TetraEdgesStiffness;
    typedef Vec<6,Mat3x3> EigenTensors;
    typedef Vec<6,Real> EigenValues;

    typedef Vec4 ParameterArray;
    // typedef Vec10 ParameterArray;

    typedef type::vector<Coord> AnisotropyDirectionArray;

    typedef core::topology::BaseMeshTopology::PointID Index;
    typedef core::topology::BaseMeshTopology::SeqTetrahedra VecElement;
    /// Rigid transformation (rotation) matrix
    typedef type::MatNoInit<3, 3, Real> Transformation;
    /// Stiffness matrix ( = RJKJtRt  with K the Material stiffness matrix, J the strain-displacement matrix, and R the transformation matrix if any )
    typedef type::Mat<12, 12, Real> StiffnessMatrix;


    /// the way the stiffness matrix should be computed on HighOrder elements
    typedef enum
    {
        ISOTROPIC=1,
        TRANSVERSE_ISOTROPIC=2,
        ORTHOTROPIC=3,
        CUBIC=4
    } ElasticitySymmetry;

protected:

    /// data structure stored for each tetrahedron
    class TetrahedronRestInformation
    {
    public:
        typedef typename DataTypes::Real  Real;
        typedef Mat<3,3,Real> Mat3x3;


        /// rest volume

        Coord shapeVector[4];
        Coord restEdgeVector[6];
        Mat3x3 rotation; // rotation from deformed to rest configuration
        Mat3x3 restRotation; // used for QR decomposition
        TetraEdgesStiffness stiffnessVector; // the stiffness matrices where is the number of control points
        EigenTensors eigenTensors; // for eigen tensors
        EigenValues eigenValues; // for eigen values
        size_t v[4]; // the indices of the 4 vertices
        TetraEdgesStiffness rotatedStiffnessVector; // the stiffness matrices where is the number of control points
        type::vector<Mat3x3> reducedStiffnessVector;
        // store 6 2x3 matrices per integration points
        type::vector<  Mat6x9 >  integrationPointsStiffnessVector;
        // store the 6 rest edge vector for each integration point
        type::vector< type::vector<Coord> > integrationPointsRestEdgeVector;
        // store a rest rotation for each integration point
        type::vector< Mat3x3 > integrationPointsRestRotationArray;

        /// Output stream
        inline friend std::ostream& operator<< ( std::ostream& os, const TetrahedronRestInformation& /*eri*/ )
        {
            return os;
        }

        /// Input stream
        inline friend std::istream& operator>> ( std::istream& in, TetrahedronRestInformation& /*eri*/ )
        {
            return in;
        }

        TetrahedronRestInformation()
        {
        }
    };

    TetrahedronData<sofa::type::vector<TetrahedronRestInformation> > tetrahedronInfo;


    sofa::core::topology::BaseMeshTopology* _topology;
    VecCoord  _initialPoints;///< the intial positions of the points

    bool updateMatrix;
    bool updateTopologyInfo;


    ElasticitySymmetry elasticitySymmetry;


    Real lambda;  /// first Lame coefficient
    Real mu;    /// second Lame coefficient



    TetrahedronAnisotropicForceField();

public:
    Data<type::vector<Real>> d_poissonRatio;
    Data<type::vector<Real>> d_youngModulus;

    Data<std::string> d_anisotropy;
    Data<type::vector<ParameterArray>> d_anisotropyParameter;
    Data<AnisotropyDirectionArray> d_anisotropyDirection;
    Data<type::vector<Vec5>> d_controlPoints;
    Data<int> d_IDWDepth;
    Data<Real> d_meshRotation;

    Data<bool> d_drawHeterogeneousTetra; ///< Draw Heterogeneous Tetra in different color
    Data<bool> d_drawDirection;
    Data<Real> d_transparency;

    void init() override;
    void reinit() override;


    virtual void addForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & dataV ) override;
    virtual void addDForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX ) override;
    virtual SReal getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const override;

//    void addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix ) override;
    void addKToMatrix(sofa::linearalgebra::BaseMatrix *mat, SReal k, unsigned int &offset) override;

    void updateTopologyInformation();
    void setMechanicalParametersFromControlPoints(size_t eltIndex,Tetra indexArray);
    type::vector<std::pair<Real, int> > generateListOfNorm(const Coord tetraBarycenter);
    void IDWdata(type::vector<std::pair<Real, int> > listOfNormWeighted,Real& dataToInterpolate,size_t indexData);
    void IDWdirection(type::vector<std::pair<Real, int> > listOfNormWeighted, Coord &dir);

    virtual Real getLambda() const { return lambda;}
    virtual Real getMu() const { return mu;}


    void draw(const core::visual::VisualParams* vparams) override;
    /// compute lambda and mu based on the Young modulus and Poisson ratio
    void updateLameCoefficients();


    void computeTetrahedronStiffnessEdgeMatrixForElts(TetraEdgesStiffness &stiffnessVector, const EigenTensors &eigenTensors, const EigenValues &eigenValues, const Coord point[4]);
    void initStiffnessVector(TetrahedronRestInformation &t, const Tetrahedron &);
    void updateStiffnessVector(const Tetrahedron &tetra, TetrahedronRestInformation &my_tinfo);

protected:

    static void computeQRRotation(Mat3x3 &r, const Coord *dp);

    virtual const TetraEdgesStiffness &getStiffnessArray(const TetrahedronRestInformation *restTetra);
    virtual const TetraEdgesStiffness &getRotatedStiffnessArray(const TetrahedronRestInformation *restTetra);

    void computeKelvinModesForElts(size_t eltIndex, EigenTensors &eigenTensors, EigenValues &eigenValues);

    /// Checks all elements for physical stability (positive Kelvin eigenvalues).
    /// Emits msg_error with per-element detail for the first unstable element and
    /// a count summary at the end. Call at the end of init() and reinit().
    void verifyStability() const;

    void sortArr(type::vector<Real> arr, int n,type::vector<std::pair<Real, int> >* vp)
    {

        // Vector to store element
        // with respective present index
        //        type::vector<std::pair<Real, int> > vp;

        // Inserting element in pair vector
        // to keep track of previous indexes
        for (int i = 0; i < n; ++i) {
            auto p = std::make_pair(arr[i], i);
            vp->push_back(p);
        }

        // Sorting pair vector
        sort(vp->begin(), vp->end());
    }

};


#if !defined(ANISOTROPY_COMPONENT_FORCEFIELD_TETRAHEDRONANISOTROPICFORCEFIELD_CPP)
extern template class SOFA_ANISOTROPY_API TetrahedronAnisotropicForceField<defaulttype::Vec3Types>;
#endif
}
