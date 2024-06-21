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
    typedef Mat<3,6,Real>       Mat3x6  ;
    typedef Mat<4,4,Real>       Mat4x4  ;
    typedef Mat<6,3,Real>       Mat6x3  ;
    typedef Mat<6,6,Real>		Mat6x6  ;
    typedef Mat<45,6,Real>		Mat45x6  ;
    typedef Mat<45,9,Real>		Mat45x9  ;
    typedef Mat<190,6,Real>     Mat190x6  ;
    typedef Mat<190,9,Real>     Mat190x9  ;
    typedef Mat<595,6,Real>     Mat595x6  ;
    typedef Mat<595,9,Real>     Mat595x9  ;
    typedef Mat<1540,6,Real>    Mat1540x6  ;
    typedef Mat<1540,9,Real>    Mat1540x9  ;
    typedef Mat<2,6,Real>       Mat2x6  ;
    typedef Mat<6,9,Real>       Mat6x9  ;


    typedef type::MatSym<3,Real> MatrixSym;
    // In case of non 3D template
    typedef Vec<3,Real> Vec3;
    typedef Vec<4,Real> Vec4;
    typedef Vec<5,Real> Vec5;
    typedef Vec<6,Real> Vec6;
    typedef Vec<9,Real> Vec9;
    typedef Vec<10,Real> Vec10;
    typedef Vec<16, Real> Vec16;
    typedef Vec<16, int> Vec16Int;
    typedef StdVectorTypes< Vec3, Vec3, Real >     GeometricalTypes ; /// assumes the geometry object type is 3D
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
           // structure that store coefficients matrices for elements of order 2,3,4,5
    // to save memory only store pointer to matrices of predefined size
    struct weightArrayPointer {
    public:
        boost::shared_ptr<Mat45x6>   weightArrayQuadratic[2];
        boost::shared_ptr<Mat190x6>  weightArrayCubic[2];
        boost::shared_ptr<Mat595x6>  weightArrayQuartic[2];
        boost::shared_ptr<Mat1540x6> weightArrayQuintic[2];

        void allocate(size_t degree);
    };
    // the array where stiffness coefficients are stored for affine elements of any degree
    std::vector<Vec6> affineStiffnessCoefficientArray;
    // the array where stiffness coefficients are stored for affine elements of degree < 5
    weightArrayPointer affineStiffnessCoefficientPreStoredArray;
    // for Bezier numerical integration store the fixed coefficients independent from the point of integration
    std::vector<Vec16> bezierCoefficientArray;
    // for Bezier numerical integration store the index mapping
    std::vector<Vec16Int> bezierMappingArray;
        // the data stored for each integration point
    struct NumericalIntegrationStiffnessData {
        // the weight of the integration point
        Real integrationWeight;
        // for each pair of control point store the weight matrix 6 w_\gamma * dN_p/d\param_i(\param_gamma) *  dN_q/d\param_j (\param_gamma)
        std::vector<Vec6> weightArray;
        // for each pair of control point store the weight matrix 6 w_\gamma * dN_p/d\param_i(\param_gamma) *  dN_q/d\param_j (\param_gamma)
        std::vector<Mat4x4> weightArray4x4;
        // for Bezier numerical integration store the coefficients that depend on the integration points
        std::vector<Real> weightBezierArray;
        weightArrayPointer arrayPointer;
        // for each control point  store the derivative of the shape functions
        std::vector<Deriv> coefficientArray;
        /// barycentric coordinate of the integration point \param_i
        Vec4 integrationPoint;
    };
    // the array where the weights and coordinates of each integration points are stored
    std::vector<NumericalIntegrationStiffnessData> numericalIntegrationStiffnessDataArray;
    /// the elasticity tensor
    Mat6x6 elasticityTensor;

    size_t nbControlPoints = 4;
    size_t nbStiffnessEntries = nbControlPoints*(nbControlPoints - 1) / 2;

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
        TetraEdgesStiffness stiffnessVector; // the nc*(nc+1)/2 stiffness matrices where nc is the number of control points
        EigenTensors eigenTensors; // for eigen tensors
        EigenValues eigenValues; // for eigen values
        size_t v[4]; // the indices of the 4 vertices
        TetraEdgesStiffness rotatedStiffnessVector; // the nc*(nc+1)/2 stiffness matrices where nc is the number of control points
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

//    class FTCFTetrahedronHandler : public TopologyDataHandler<Tetrahedron, sofa::type::vector<TetrahedronRestInformation> >
//    {
//    public:
//        typedef typename TetrahedronAnisotropicForceField<DataTypes>::TetrahedronRestInformation TetrahedronRestInformation;

//        FTCFTetrahedronHandler(TetrahedronAnisotropicForceField<DataTypes>* ff,
//                               TetrahedronData<sofa::type::vector<TetrahedronRestInformation> >* data )
//            :TopologyDataHandler<Tetrahedron, sofa::type::vector<TetrahedronRestInformation> >(data)
//            ,ff(ff)
//        {

//        }

//        void applyCreateFunction(unsigned int, TetrahedronRestInformation &t, const Tetrahedron
//                                                                                  &, const sofa::type::vector<unsigned int> &, const sofa::type::vector<double> &);

//        void updateStiffnessVector(unsigned int tetrahedronIndex,TetrahedronRestInformation &my_tinfo);
////        void interpolateBewteenControlPoints(unsigned int interpolateBewteenControlPoints,TetrahedronRestInformation &my_tinfo);
//    protected:
//        TetrahedronAnisotropicForceField<DataTypes>* ff;

//    };

    TetrahedronData<sofa::type::vector<TetrahedronRestInformation> > tetrahedronInfo;


    sofa::core::topology::BaseMeshTopology* _topology;
    VecCoord  _initialPoints;///< the intial positions of the points

    bool updateMatrix;
    bool updateTopologyInfo;


    ElasticitySymmetry elasticitySymmetry;


    Real lambda;  /// first Lame coefficient
    Real mu;    /// second Lame coefficient



    TetrahedronAnisotropicForceField();

//    virtual ~TetrahedronAnisotropicForceField();

public:
    // user input
    Data<std::string> d_method; ///< the computation method of the displacements


    Data<type::vector<Real>> d_poissonRatio; // stiffness coefficient for isotropic elasticity;
    Data<type::vector<Real>> d_youngModulus;

    Data<std::string> d_anisotropy; // the type of isotropy
    Data<type::vector<ParameterArray>> d_anisotropyParameter; // the set of parameters defining the elasticity anisotropy
    Data<AnisotropyDirectionArray> d_anisotropyDirection; // the directions of anisotropy
//    Data<type::vector<type::vector<Mat3x3>>> d_ortho_matrix;
//    Data<type::vector<type::vector<Real>>> d_ortho_scalar;
    Data<type::vector<Vec5>> d_controlPoints;
    Data<int> d_IDWDepth;
    Data<Real> d_meshRotation;
    //    Data<Real> d_poissonRatio; // stiffness coefficient for isotropic elasticity;
    //    Data<Real> d_youngModulus;

    //	Data<std::string> d_anisotropy; // the type of isotropy
    //    Data<ParameterArray> d_anisotropyParameter; // the set of parameters defining the elasticity anisotropy
    //	Data<AnisotropyDirectionArray> d_anisotropyDirection; // the directions of anisotropy

    /// the order of integration for numerical integration
//    Data<size_t>	     numericalIntegrationOrder;
    /// the type of numerical integration method chosen
//    Data<std::string>	     numericalIntegrationMethod;
    /// the type of integration method chosen for non linear element.
//    Data<std::string>	 d_integrationMethod;
    // if one rotation is attached to each integration point
    Data<bool> d_oneRotationPerIntegrationPoint;
    // measure the time spent in assembling the stiffness matrix
//    Data<Real> d_assemblyTime;
    // whether each affine element should be assembled with the affine assembly method irrespective to the chosen integration method
//    Data<bool> d_forceAffineAssemblyForAffineElements;

    Data< bool > d_drawHeterogeneousTetra; ///< Draw Heterogeneous Tetra in different color
    Data< bool > d_drawDirection;
    Data< Real > d_transparency;
//    Data< bool > test_visu;

    Real m_yaw;
    Real m_roll;

    void init() override;
    void reinit() override;


    virtual void addForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & dataV ) override;
    virtual void addDForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX ) override;
    virtual SReal getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const override;

//    void addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix ) override;
    void addKToMatrix(sofa::linearalgebra::BaseMatrix *mat, SReal k, unsigned int &offset) override;

    void updateTopologyInformation();
    void setMechanicalParametersFromControlPoints(size_t eltIndex,Tetra indexArray);
    type::vector<std::pair<Real, int> > generateListOfNorm(Coord tetraBarycenter);
    void IDWInterpolationBewteenControlPoints(type::vector<std::pair<Real, int> > listOfNormWeighted,Real* dataToInterpolate,size_t indexData);

    virtual Real getLambda() const { return lambda;}
    virtual Real getMu() const { return mu;}


    void draw(const core::visual::VisualParams* vparams) override;
    /// compute lambda and mu based on the Young modulus and Poisson ratio
    void updateLameCoefficients();


//    void computeTetrahedronStiffnessEdgeMatrix(const Coord position[4],Mat6x9 edgeStiffness[2]);
//    void computeTetrahedronStiffnessEdgeMatrix(const Coord position[4],Mat3x3 edgeStiffness[6]);

    /////////////// TEST ////////////////////////
    void computeTetrahedronStiffnessEdgeMatrixForElts(TetraEdgesStiffness &stiffnessVector,const EigenTensors & eigenTensors , const EigenValues & eigenValues, const Coord point[4]);
//    void updateStiffnessVectorWithCP();
//    void IDWInterpolationBewteenControlPoints();
    /////////////////////////////////////////////

    void initStiffnessVector(TetrahedronRestInformation &t, const Tetrahedron &);

    void updateStiffnessVector(const Tetrahedron &tetra,TetrahedronRestInformation &my_tinfo);
//        void interpolateBewteenControlPoints(unsigned int interpolateBewteenControlPoints,TetrahedronRestInformation &my_tinfo);

//    friend class FTCFTetrahedronHandler;

protected :

//    FTCFTetrahedronHandler* tetrahedronHandler;

    static void computeQRRotation( Mat3x3 &r, const Coord *dp);

    virtual const TetraEdgesStiffness &getStiffnessArray(const TetrahedronRestInformation *restTetra);

    virtual const TetraEdgesStiffness &getRotatedStiffnessArray(const TetrahedronRestInformation *restTetra);

    // compute elasticity tensor for isotropic and anisotropic cases
//    void computeElasticityTensor();

//    void computeKelvinModes();
    void computeKelvinModesForElts(size_t eltIndex,EigenTensors &eigenTensors , EigenValues &eigenValues);
    //    Mat3x3 matrixD(Real _lambda,Vec3 v1,Vec3 v2,Vec3 n,Eigen::MatrixXd m);
    //    Mat3x3 matrixZ(Real _lambda_k, Real _lambda_i, Real _lambda_j,Eigen::MatrixXd m);

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
