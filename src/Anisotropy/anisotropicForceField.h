
#ifndef SOFA_COMPONENT_FORCEFIELD_HIGHORDERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_HIGHORDERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif
#include "initHighOrderFEM.h"
#include <sofa/core/behavior/ForceField.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/helper/fixed_array.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <SofaBaseTopology/TopologyData.h>
#include <sofa/defaulttype/MatSym.h>
#include <sofa/helper/system/thread/CTime.h>
#include <sofa/defaulttype/SolidTypes.h>

namespace sofa
{

namespace component
{
namespace topology
{
	template< class DataTypes> class HighOrderTetrahedronSetGeometryAlgorithms;
}
namespace forcefield
{

using namespace sofa::defaulttype;
using namespace sofa::component::topology;



template<class DataTypes>
class SOFA_HIGHORDER_FEM_API HighOrderTetrahedralCorotationalFEMForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HighOrderTetrahedralCorotationalFEMForceField,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::Real        Real        ;
    typedef typename DataTypes::Coord       Coord       ;
    typedef typename DataTypes::Deriv       Deriv       ;
    typedef typename DataTypes::VecCoord    VecCoord    ;
    typedef typename DataTypes::VecDeriv    VecDeriv    ;
    typedef typename DataTypes::VecReal     VecReal     ;
    typedef Data<VecCoord>                  DataVecCoord;
    typedef Data<VecDeriv>                  DataVecDeriv;    
	typedef sofa::core::topology::Topology::Tetrahedron Tetrahedron;
	typedef sofa::core::topology::Topology::TetraID TetraID;
	typedef sofa::core::topology::Topology::Tetra Tetra;
	typedef sofa::core::topology::Topology::Point Point;
	typedef sofa::core::topology::Topology::Triangle Triangle;
	typedef sofa::core::topology::Topology::Edge Edge;
	typedef sofa::core::topology::Topology::Quad Quad;
//        typedef helper::Quater<Real> Rot;
	typedef sofa::core::topology::Topology::Hexahedron Hexahedron;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInTriangle EdgesInTriangle;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInTetrahedron EdgesInTetrahedron;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInQuad EdgesInQuad;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInHexahedron EdgesInHexahedron;
	typedef sofa::core::topology::BaseMeshTopology::TrianglesInTetrahedron TrianglesInTetrahedron;
	typedef sofa::core::topology::BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;

	typedef helper::vector<Coord> SetAnisotropyDirectionArray; // When the model is anisotropic, for instance in invariant I4
	typedef helper::vector<Real> SetParameterArray; //necessary to store hyperelastic parameters (mu, lambda ...)

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
//	typedef Mat<9,6,Real>       Mat9x6  ;
	typedef Mat<6,9,Real>       Mat6x9  ;


	typedef typename defaulttype::MatSym<3,Real> MatrixSym;
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
    typedef typename sofa::component::topology::HighOrderTetrahedronSetGeometryAlgorithms<GeometricalTypes>::VecPointID VecPointID;
    typedef typename defaulttype::SolidTypes<Real>::Transform       Transform;

    typedef Vec4 ParameterArray;
    // typedef Vec10 ParameterArray;
    
    typedef helper::vector<Coord> AnisotropyDirectionArray;

    typedef core::topology::BaseMeshTopology::PointID Index;
    typedef core::topology::BaseMeshTopology::SeqTetrahedra VecElement;
    /// Rigid transformation (rotation) matrix
    typedef defaulttype::MatNoInit<3, 3, Real> Transformation;
    /// Stiffness matrix ( = RJKJtRt  with K the Material stiffness matrix, J the strain-displacement matrix, and R the transformation matrix if any )
    typedef defaulttype::Mat<12, 12, Real> StiffnessMatrix;

    typedef enum
    {
        POLAR_DECOMPOSITION,
        QR_DECOMPOSITION,
	    POLAR_DECOMPOSITION_MODIFIED,
		LINEAR_ELASTIC
    } RotationDecompositionMethod;

	/// the way the stiffness matrix should be computed on HighOrder elements
    typedef enum
    {
        AFFINE_ELEMENT_INTEGRATION = 1,
        NUMERICAL_INTEGRATION = 2,
        STANDARD_INTEGRATION = 3,
        BEZIER_NUMERICAL_INTEGRATION = 4,
        NUMERICAL_INTEGRATION_2 = 5
    } IntegrationMethod;

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
		helper::vector<Mat3x3> stiffnessVector; // the nc*(nc+1)/2 stiffness matrices where nc is the number of control points
		size_t v[4]; // the indices of the 4 vertices
		helper::vector<Mat3x3> rotatedStiffnessVector; // the nc*(nc+1)/2 stiffness matrices where nc is the number of control points
		helper::vector<Mat3x3> reducedStiffnessVector;
		// store 6 2x3 matrices per integration points
		helper::vector<  Mat6x9 >  integrationPointsStiffnessVector;
		// store the 6 rest edge vector for each integration point
		helper::vector< helper::vector<Coord> > integrationPointsRestEdgeVector;
		// store a rest rotation for each integration point
		helper::vector< Mat3x3 > integrationPointsRestRotationArray;

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

    class FTCFTetrahedronHandler : public TopologyDataHandler<Tetrahedron, sofa::helper::vector<TetrahedronRestInformation> >
    {
    public:
        typedef typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation TetrahedronRestInformation;

        FTCFTetrahedronHandler(HighOrderTetrahedralCorotationalFEMForceField<DataTypes>* ff,
                TetrahedronData<sofa::helper::vector<TetrahedronRestInformation> >* data )
            :TopologyDataHandler<Tetrahedron, sofa::helper::vector<TetrahedronRestInformation> >(data)
            ,ff(ff)
        {

        }

        void applyCreateFunction(unsigned int, TetrahedronRestInformation &t, const Tetrahedron
                &, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &);

        void updateStiffnessVector(unsigned int tetrahedronIndex,TetrahedronRestInformation &my_tinfo);
        void interpolateBewteenControlPoints(unsigned int interpolateBewteenControlPoints,TetrahedronRestInformation &my_tinfo);
    protected:
        HighOrderTetrahedralCorotationalFEMForceField<DataTypes>* ff;

    };

    TetrahedronData<sofa::helper::vector<TetrahedronRestInformation> > tetrahedronInfo;


    sofa::core::topology::BaseMeshTopology* _topology;
	sofa::component::topology::HighOrderTetrahedronSetGeometryAlgorithms<GeometricalTypes>* highOrderTetraGeo;
    VecCoord  _initialPoints;///< the intial positions of the points

    bool updateMatrix;
    bool updateTopologyInfo;


    RotationDecompositionMethod decompositionMethod;

	ElasticitySymmetry elasticitySymmetry;

	std::vector<Mat3x3> anisotropyMatrixArray;
	std::vector<Real> anisotropyScalarArray;

    /////////////// TEST ////////////////////////
    std::vector<std::vector<Mat3x3>> vecAnisotropyMatrixArray;
    std::vector<std::vector<Real>> vecAnisotropyScalarArray;
    ////////////////////////////////////////////

    Real lambda;  /// first Lame coefficient
    Real mu;    /// second Lame coefficient

	
	IntegrationMethod    integrationMethod;



    HighOrderTetrahedralCorotationalFEMForceField();

    virtual ~HighOrderTetrahedralCorotationalFEMForceField();

public:
	// user input
	Data<std::string> d_method; ///< the computation method of the displacements


    Data<helper::vector<Real>> d_poissonRatio; // stiffness coefficient for isotropic elasticity;
    Data<helper::vector<Real>> d_youngModulus;

    Data<std::string> d_anisotropy; // the type of isotropy
    Data<helper::vector<ParameterArray>> d_anisotropyParameter; // the set of parameters defining the elasticity anisotropy
    Data<AnisotropyDirectionArray> d_anisotropyDirection; // the directions of anisotropy
    Data<helper::vector<helper::vector<Mat3x3>>> d_ortho_matrix;
    Data<helper::vector<helper::vector<Real>>> d_ortho_scalar;
    Data<helper::vector<Vec5>> d_controlPoints;
    Data<Real> d_meshRotation;
//    Data<Real> d_poissonRatio; // stiffness coefficient for isotropic elasticity;
//    Data<Real> d_youngModulus;

//	Data<std::string> d_anisotropy; // the type of isotropy
//    Data<ParameterArray> d_anisotropyParameter; // the set of parameters defining the elasticity anisotropy
//	Data<AnisotropyDirectionArray> d_anisotropyDirection; // the directions of anisotropy

	/// the order of integration for numerical integration
	Data<size_t>	     numericalIntegrationOrder;
	/// the type of numerical integration method chosen
	Data<std::string>	     numericalIntegrationMethod;
	/// the type of integration method chosen for non linear element.
	Data<std::string>	 d_integrationMethod; 
	// if one rotation is attached to each integration point 
	Data<bool> d_oneRotationPerIntegrationPoint;
	// measure the time spent in assembling the stiffness matrix
	Data<Real> d_assemblyTime;
	// whether each affine element should be assembled with the affine assembly method irrespective to the chosen integration method  
	Data<bool> d_forceAffineAssemblyForAffineElements;

    Data< bool > d_drawHeterogeneousTetra; ///< Draw Heterogeneous Tetra in different color
    Data< bool > d_drawDirection;
    Data< Real > d_transparency;
    Data< bool > test_visu;

    Real m_yaw;
    Real m_roll;

    void init() override;
    void reinit() override;


    virtual void addForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & dataV ) ;
    virtual void addDForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX ) ;
    virtual SReal getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const;

    void addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix ) override;
    void addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal k, unsigned int &offset) override;

    void updateTopologyInformation();

    virtual Real getLambda() const { return lambda;}
    virtual Real getMu() const { return mu;}


    void setRotationDecompositionMethod( const RotationDecompositionMethod m)
    {
        decompositionMethod=m;
    }
    void draw(const core::visual::VisualParams* vparams);
    /// compute lambda and mu based on the Young modulus and Poisson ratio
    void updateLameCoefficients();


    void computeTetrahedronStiffnessEdgeMatrix(const Coord position[4],Mat6x9 edgeStiffness[2]);
    void computeTetrahedronStiffnessEdgeMatrix(const Coord position[4],Mat3x3 edgeStiffness[6]);

    /////////////// TEST ////////////////////////
    void computeTetrahedronStiffnessEdgeMatrixForElts(size_t eltIndex,const Coord point[4], Mat6x9 edgeStiffnessVectorized[2]);
    void updateStiffnessVectorWithCP();
    void IDWInterpolationBewteenControlPoints();
    /////////////////////////////////////////////

	friend class FTCFTetrahedronHandler;

protected :
    FTCFTetrahedronHandler* tetrahedronHandler;
	// for profiling the assembly task
	helper::system::thread::ctime_t totalUpdateMat;
	helper::system::thread::ctime_t totalComputeLocalStiffness;
	static void computeQRRotation( Mat3x3 &r, const Coord *dp);

	virtual const helper::vector<Mat3x3> &getStiffnessArray(const size_t i,
		TetrahedronRestInformation *restTetra);

	virtual const helper::vector<Mat3x3> &getRotatedStiffnessArray(const size_t i,
		const TetrahedronRestInformation *restTetra);

	// compute elasticity tensor for isotropic and anisotropic cases
    void computeElasticityTensor();

    void computeKelvinModes();
    void computeKelvinModesForElts(size_t eltIndex);
//    Mat3x3 matrixD(Real _lambda,Vec3 v1,Vec3 v2,Vec3 n,Eigen::MatrixXd m);
//    Mat3x3 matrixZ(Real _lambda_k, Real _lambda_i, Real _lambda_j,Eigen::MatrixXd m);

    void sortArr(helper::vector<Real> arr, int n,helper::vector<std::pair<Real, int> >* vp)
    {

        // Vector to store element
        // with respective present index
//        helper::vector<std::pair<Real, int> > vp;

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


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_FASTTETRAHEDRALCOROTATIONALFORCEFIELD_H