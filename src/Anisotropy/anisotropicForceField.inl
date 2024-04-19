
template< class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::FTCFTetrahedronHandler::applyCreateFunction(unsigned int tetrahedronIndex,
    TetrahedronRestInformation &my_tinfo,
    const Tetrahedron &tt,
    const sofa::helper::vector<unsigned int> &,
    const sofa::helper::vector<double> &)
{
    if (ff)
    {


        const std::vector< Tetrahedron > &tetrahedronArray = ff->_topology->getTetrahedra();
        HighOrderTetrahedronSetTopologyContainer *container = ff->highOrderTetraGeo->getTopologyContainer();
        HighOrderDegreeType degree = container->getDegree();
        size_t nbControlPoints = (degree + 1)*(degree + 2)*(degree + 3) / 6;
        size_t nbStiffnessEntries = nbControlPoints*(nbControlPoints - 1) / 2;
        if (my_tinfo.stiffnessVector.size() != nbStiffnessEntries) {
            my_tinfo.stiffnessVector.resize(nbStiffnessEntries);
        }
        // set array to zero
        std::fill(my_tinfo.stiffnessVector.begin(), my_tinfo.stiffnessVector.end(), Mat3x3());

        //      const std::vector< Edge> &edgeArray=ff->_topology->getEdges() ;
        size_t i, j, k, l, m, n;

        typename DataTypes::Coord point[4];


        const typename DataTypes::VecCoord &restPosition = ff->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
        // now computed the stiffness for the HighOrder Tetrahedron
        sofa::helper::vector<TetrahedronIndexVector> tbiArray;

        tbiArray = ff->highOrderTetraGeo->getTopologyContainer()->getTetrahedronIndexArray();

        size_t rank = 0;


        ///describe the indices of the 4 tetrahedron vertices
        const Tetrahedron &t = tetrahedronArray[tetrahedronIndex];
        //    BaseMeshTopology::EdgesInTetrahedron te=ff->_topology->getEdgesInTetrahedron(tetrahedronIndex);


        // store the point position
        for (j = 0; j < 4; ++j)
            point[j] = (restPosition)[t[j]];

        if (ff->decompositionMethod == QR_DECOMPOSITION) {
            for (j = 0; j < 6; ++j) {
                k = edgesInTetrahedronArray[j][0];
                l = edgesInTetrahedronArray[j][1];

                // store the rest edge vector
                my_tinfo.restEdgeVector[j] = point[l] - point[k];
            }
            // compute the rotation matrix of the initial tetrahedron for the QR decomposition
            computeQRRotation(my_tinfo.restRotation, my_tinfo.restEdgeVector);
        }
        else    if (ff->decompositionMethod == POLAR_DECOMPOSITION_MODIFIED) {
            Mat3x3 Transformation;
            Transformation[0] = point[1] - point[0];
            Transformation[1] = point[2] - point[0];
            Transformation[2] = point[3] - point[0];
            helper::Decompose<Real>::polarDecomposition(Transformation, my_tinfo.restRotation);
        }

        if ((ff->integrationMethod == HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::AFFINE_ELEMENT_INTEGRATION) ||
            ((ff->d_forceAffineAssemblyForAffineElements.getValue()) && ((ff->highOrderTetraGeo->isBezierTetrahedronAffine(tetrahedronIndex, restPosition))))) {
            Mat6x9 edgeStiffnessVectorized[2];

            if (ff->d_anisotropyDirection.getValue().size() == ff->_topology->getNbTetrahedra())
                ff->computeTetrahedronStiffnessEdgeMatrixForElts(tetrahedronIndex,point, edgeStiffnessVectorized);
            else
                ff->computeTetrahedronStiffnessEdgeMatrix(point, edgeStiffnessVectorized);

            helper::system::thread::ctime_t startUpdateMat = helper::system::thread::CTime::getTime();
            if (degree == 1) {
                for (rank = 0; rank < nbStiffnessEntries; rank++) {
                    my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&edgeStiffnessVectorized[0][rank][0]);
                }
            }

            }

            if (ff->f_printLog.getValue()) {
                helper::system::thread::ctime_t endUpdateMat = helper::system::thread::CTime::getTime();
                ff->totalUpdateMat += endUpdateMat - startUpdateMat;
            }

        }
    }


template <class DataTypes>
HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::HighOrderTetrahedralCorotationalFEMForceField()
        : tetrahedronInfo(initData(&tetrahedronInfo, "tetrahedronInfo", "Internal tetrahedron data"))
        , _initialPoints(0)
        , updateMatrix(true)
        , d_method(initData(&d_method, std::string("linear"), "method", "method for rotation computation :\"qr\" (by QR) or \"polar\" or \"polar2\" or \"none\" (Linear elastic)"))

        , d_poissonRatio(initData(&d_poissonRatio, "poissonRatio", "Poisson ratio in Hooke's law"))
        , d_youngModulus(initData(&d_youngModulus, "youngModulus", "Young modulus in Hooke's law"))
        , d_anisotropy(initData(&d_anisotropy, std::string("isotropic"), "elasticitySymmetry", "the type of anisotropy for the elasticity tensor :\"isotropic\"  or \"transverseIsotropic\" or \"ortho_scalar\" or \"cubic\" "))
        , d_anisotropyParameter(initData(&d_anisotropyParameter, "anisotropyParameters", "the elastic parameters for anisotropic materials.\n"
                                                                                         "- for cubic symmetry       --> anisotropyParameters == [anisotropyRatio]\n"
                                                                                         "- for transverse symemetry --> anisotropyParameters == [youngModulusLongitudinal, poissonRatioTransverseLongitudinal, shearModulusTransverse]"))
        , d_ortho_scalar(initData(&d_ortho_scalar,"ortho_scalar",""))
        , d_ortho_matrix(initData(&d_ortho_matrix,"ortho_matrix",""))
        , d_anisotropyDirection(initData(&d_anisotropyDirection, "anisotropyDirections", "the directions of anisotropy"))
        , numericalIntegrationOrder(initData(&numericalIntegrationOrder, (size_t)2, "integrationOrder", "The order of integration for numerical integration"))
        , d_integrationMethod(initData(&d_integrationMethod, std::string("analytical"), "integrationMethod", "\"analytical\" if closed form expression for affine element, \"numerical\" if numerical integration is chosen,  \"standard\" if standard integration is chosen"))
        , numericalIntegrationMethod(initData(&numericalIntegrationMethod, std::string("Tetrahedron Gauss"), "numericalIntegrationMethod", "The type of numerical integration method chosen"))
        , d_oneRotationPerIntegrationPoint(initData(&d_oneRotationPerIntegrationPoint, false, "oneRotationPerIntegrationPoint", "if true then computes one rotation per integration point"))
        , d_assemblyTime(initData(&d_assemblyTime, (Real)0, "assemblyTime", "the time spent in assembling the stiffness matrix. Only updated if printLog is set to true"))
        , d_forceAffineAssemblyForAffineElements(initData(&d_forceAffineAssemblyForAffineElements, true, "forceAffineAssemblyForAffineElements", "if true affine tetrahedra are always assembled with the closed form formula, Otherwise use the method defined in integrationMethod"))
        , d_drawHeterogeneousTetra(initData(&d_drawHeterogeneousTetra,true,"drawHeterogeneousTetra","Draw Heterogeneous Tetra in different color"))
        , d_drawDirection(initData(&d_drawDirection,true,"drawDirection","Draw different color for each direction"))
        , d_transparency(initData(&d_transparency,(Real)0.25,"transparency","transparency [0,1]"))
        , d_controlPoints(initData(&d_controlPoints,"controlPoints","controlPoints"))
        , d_meshRotation(initData(&d_meshRotation,"meshRotation",""))
        , m_yaw(45)
        , m_roll(45)
        , lambda(0)
        , mu(0)
        , tetrahedronHandler(NULL)
    {
        tetrahedronHandler = new FTCFTetrahedronHandler(this, &tetrahedronInfo);
    }


template <class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::init()
{
    this->d_componentState = ComponentState::Invalid;

    //  serr << "initializing HighOrderTetrahedralCorotationalFEMForceField" << sendl;
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();
    this->getContext()->get(highOrderTetraGeo);


    if (d_anisotropy.getValue() == "transverseIsotropic") 
        elasticitySymmetry= TRANSVERSE_ISOTROPIC;


    if ((d_method.getValue() == "qr") || (d_method.getValue() == "large"))
        decompositionMethod= QR_DECOMPOSITION;

    if (d_integrationMethod.getValue() == "analytical")
        integrationMethod = AFFINE_ELEMENT_INTEGRATION;


    helper::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    tetrahedronInf.resize(_topology->getNbTetrahedra());

    if (_initialPoints.size() == 0)
    {
        // get restPosition
        const VecCoord& p = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
        _initialPoints=p;
    }

    helper::system::thread::ctime_t startAssembly=helper::system::thread::CTime::getTime();
    totalUpdateMat=0;
    totalComputeLocalStiffness=0;

    size_t i;


    // precompute the coefficients for handling affine elements
    topology::TetrahedronIndexVector tbi1,tbi2;
    affineStiffnessCoefficientArray.clear();
    
    topology::HighOrderDegreeType degree=highOrderTetraGeo->getTopologyContainer()->getDegree();
    size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
    sofa::helper::vector<topology::TetrahedronIndexVector> tbiArray;
    tbiArray=highOrderTetraGeo->getTopologyContainer()->getTetrahedronIndexArray();

    if (degree==1)
    {
        msg_warning() << "The order of your TopologyContainer is 1. With this order we force the integrationMethod to analytical";
        integrationMethod= AFFINE_ELEMENT_INTEGRATION;

        if (d_controlPoints.getValue().size() > 0)
        {
            msg_info() << "TEST INIT WITH CONTROL POINTS";

            // init array
            vecAnisotropyMatrixArray.resize(_topology->getNbTetrahedra());
            vecAnisotropyScalarArray.resize(_topology->getNbTetrahedra());
//            interpolateBewteenControlPoints();
            IDWInterpolationBewteenControlPoints();

            for (i=0; i<_topology->getNbTetrahedra(); ++i)
            {
                computeKelvinModesForElts(i);
            }
        }
        else if (d_anisotropyDirection.getValue().size() == _topology->getNbTetrahedra())
        {
            msg_info() << "TEST PUT ANISOTROPY DIRECTION TO EACH ELT OF THE MESH";
            // init array
            vecAnisotropyMatrixArray.resize(_topology->getNbTetrahedra());
            vecAnisotropyScalarArray.resize(_topology->getNbTetrahedra());

            for (i=0; i<_topology->getNbTetrahedra(); ++i)
            {
                computeKelvinModesForElts(i);
            }
        }
        else
            computeKelvinModes();
    }


    for (i=0; i<_topology->getNbTetrahedra(); ++i)
    {
        tetrahedronHandler->applyCreateFunction(i,tetrahedronInf[i],_topology->getTetrahedron(i),
                (const helper::vector< unsigned int > )0,
                (const helper::vector< double >)0);
    }

    //msg_info() << getStiffnessArray(3,&tetrahedronInf[3]);

    helper::system::thread::ctime_t endComputeLocalStiffness=helper::system::thread::CTime::getTime();
    if (this->f_printLog.getValue()){
        helper::system::thread::ctime_t endAssembly=helper::system::thread::CTime::getTime();
        std::cerr<< "Assembly time="<< ((endAssembly-startAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
        std::cerr<<" total update mat="<<((totalUpdateMat)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
        
        std::cerr<<" total compute local stiffness ="<<((totalComputeLocalStiffness)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
        
        d_assemblyTime.setValue(((endAssembly-startAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()));

        //startAssembly=helper::system::thread::CTime::getTime();
        //helper::system::thread::CTime::sleep(1.0);
        //endAssembly=helper::system::thread::CTime::getTime();
        //std::cerr<< "test time="<< ((endAssembly-startAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
    }
    /// set the call back function upon creation of a tetrahedron
    tetrahedronInfo.createTopologyHandler(_topology,tetrahedronHandler);
    tetrahedronInfo.registerTopologicalData();
    tetrahedronInfo.endEdit();

    updateTopologyInfo=true;
    this->d_componentState = ComponentState::Valid ;

}


Mat<3,3,double> matrixD(double _lambda,Vec<3,double> v1,Vec<3,double> v2,Vec<3,double> n,Eigen::Matrix3d m)
{
    double val1_D = m(0,1)*m(1,2) - m(0,1) * (m(1,1) - _lambda);
    double val2_D = m(0,1)*m(0,2) - m(1,2) * (m(0,0) - _lambda);
    double val3_D = (m(0,0) - _lambda)*(m(1,1) - _lambda) - pow(m(0,1),2);

    return val1_D * defaulttype::dyad(v1,v1) + val2_D * defaulttype::dyad(v2,v2) + val3_D * defaulttype::dyad(n,n);
}

double matrixZ(double _lambda_k, double _lambda_i, double _lambda_j,Eigen::Matrix3d m)
{
    return 1 / ( (_lambda_k - _lambda_i)*(_lambda_k - _lambda_j)*(m(0,1)*(pow(m(0,2),2) - pow(m(1,2),2)) - m(0,2)*m(1,2)*(m(0,0) - m(1,1)))
                 * ( (m(0,1)*m(0,2) - m(1,2)*(m(0,0) - _lambda_i))*( (m(0,0) - _lambda_j)  + m(0,1) + m(0,2)) )
                 - (m(0,1)*m(1,2) - m(0,2)*(m(1,1) - _lambda_i))*(m(0,1) + (m(1,1) - _lambda_j) + m(1,2)) );
}


template<class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::computeKelvinModesForElts(size_t eltIndex)
{
//    msg_info() << "-------------------------------- ENTER computeKelvinModesForElts";

    if (elasticitySymmetry != ISOTROPIC)
    {
        Vec4 anisotropyParameter=d_anisotropyParameter.getValue()[eltIndex];
        // Vec10 anisotropyParameter=d_anisotropyParameter.getValue()[eltIndex];
        //msg_info() << elasticitySymmetry;

        Real youngModulus = d_youngModulus.getValue()[eltIndex];
        Real poissonRatio = d_poissonRatio.getValue()[eltIndex];
//        if (eltIndex == 0) msg_info() << youngModulus << "  " <<anisotropyParameter;

        Coord n = d_anisotropyDirection.getValue()[eltIndex];

        n/=n.norm();
        Coord v1,v2;
        if ((n[0]!=0) || (n[1]!=0)) {
            v1=Coord(-n[1],n[0],n[2]);
        } else {
            v1=Coord(1,0,0);
        }
        v1=cross(n,v1);
        v1/=v1.norm();
        v2=cross(v1,n);

        // Here we will use the work of Sandrine Germian see thesis (2015): https://opus4.kobv.de/opus4-fau/frontdoor/index/index/docId/3490
        // Chap 4 Spectral decomposition and the Kelvin modes :
        //
        // - For CUBIC                  -------> 4.4.1 The cubic crystal system Materials (p38)
        //
        // - For TRANSVERSE_ISOTROPIC   -------> 4.4.6 The tetragonal crystal system Materials (p48)
        //
        // As well as the equations from http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm

        if (anisotropyParameter[0]==CUBIC){
            // get the different constants : anisotropy ratio.
            Real anisotropyRatio=anisotropyParameter[1];

            // Same as the isotropic but with a coefficient of anysotropy (we will call A)
            // if A == 1    --> the material is isotropic
            // defined by the expression :
            //
            // mu = A * E/(2*(1+v))
            //
            //
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //                        | / | 0   | 1   | 2   | 3        | 4        | 5         |
            //                        +===+=====+=====+=====+==========+==========+===========+
            //                        | 0 | 1-v | v   | v   | 0        | 0        | 0         |
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //                        | 1 | v   | 1-v | v   | 0        | 0        | 0         |
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //  E/(1+v)(1-2v) *       | 2 | v   | v   | 1-v | 0        | 0        | 0         |
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //                        | 3 | 0   | 0   | 0   | A(1-2v)/2 | 0        | 0        |
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //                        | 4 | 0   | 0   | 0   | 0        | A(1-2v)/2 | 0        |
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //                        | 5 | 0   | 0   | 0   | 0        | 0        | A(1-2v)/2 |
            //                        +---+-----+-----+-----+----------+----------+-----------+

            Real c11 = youngModulus * (1 - poissonRatio) / (1 - poissonRatio - 2 * poissonRatio * poissonRatio);
            Real c12 = youngModulus * poissonRatio / (1 - poissonRatio - 2 * poissonRatio * poissonRatio);
            Real c44 = anisotropyRatio * (c11 - c12) / 2;

            // The number of modes is equal to three, so that the tensor has three eigenvalues:
            Real eigen1 = c11 + 2 * c12;    // (4.39) dim 1
            Real eigen2 = c11 - c12;        // (4.40) dim 2
            Real eigen3 = 2 * c44;          // (4.41) dim 3

            // ----------------------------------------------------------------------
            // BUILD THE ORTHOGONAL MATRICES
            // ----------------------------------------------------------------------

            // The spectral decomposition of the cubic tensor consists of a dilatation,
            // a two-dimensional and a three-dimensional eigenspace
            //  (here 'x' is the dyadic product)

            // The first projection tensor associated with the first eigenvalue depends on the dilatation mode as for the isotropic tensor
            // P1 = Nd x Nd

            Mat3x3 Nd;
            Nd.identity();
            Nd/= sqrt(3.0f);

            // The second projection tensor associated with the second eigenvalue is the sum of the tensor product of the isochoric extension and pure shear modes with themselves
            // P2 = Nei x Nei + Npi x Npi (with i == 1/2/3)

            Mat3x3 Ne=(2*defaulttype::dyad(n,n)-defaulttype::dyad(v1,v1)-defaulttype::dyad(v2,v2))/sqrt(6.0f);
            Mat3x3 Np=(defaulttype::dyad(v1,v1)-defaulttype::dyad(v2,v2))/sqrt(2.0f);

            // The third projection tensor associated with the third eigenvalue is the sum of the tensor product of the three simple shear modes with themselves
            // P3 = Ns1 x Ns1 + Ns2 x Ns2 + Ns3 x Ns3

            Mat3x3 Ns1=(defaulttype::dyad(v2,n)+defaulttype::dyad(n,v2))/sqrt(2.0f);
            Mat3x3 Ns2=(defaulttype::dyad(v1,n)+defaulttype::dyad(n,v1))/sqrt(2.0f);
            Mat3x3 Ns3=(defaulttype::dyad(v1,v2)+defaulttype::dyad(v2,v1))/sqrt(2.0f);


            // push all symmetric matrices and the eigenvalues
            std::vector<Mat3x3> anisotropyMatrixArray;
            std::vector<Real> anisotropyScalarArray;

            anisotropyMatrixArray.push_back(Nd);
            anisotropyScalarArray.push_back(eigen1);
            anisotropyMatrixArray.push_back(Ne);
            anisotropyScalarArray.push_back(eigen2);
            anisotropyMatrixArray.push_back(Np);
            anisotropyScalarArray.push_back(eigen2);
            anisotropyMatrixArray.push_back(Ns1);
            anisotropyScalarArray.push_back(eigen3);
            anisotropyMatrixArray.push_back(Ns2);
            anisotropyScalarArray.push_back(eigen3);
            anisotropyMatrixArray.push_back(Ns3);
            anisotropyScalarArray.push_back(eigen3);

            vecAnisotropyMatrixArray.push_back(anisotropyMatrixArray);
            vecAnisotropyScalarArray.push_back(anisotropyScalarArray);

            //msg_info() << eigen1 << ' ' << eigen2 << ' '<< eigen2 << ' '<< eigen3 << ' ' << eigen3 << ' '<< eigen3 << ' ';
            //msg_info() << Nd << ' ' << Ne<< ' '<< Np<< ' '<< Ns1<< ' ' << Ns2<< ' '<< Ns3<< ' ';
            //msg_info() << "Nd " << Nd;
            //msg_info() << "Ne " << Ne;
            //msg_info() << "Np " << Np;
            //msg_info() << "Ns1 " << Ns1;
            //msg_info() << "Ns2 " << Ns2;
            //msg_info() << "Ns3 " << Ns3;
        }
        else if (anisotropyParameter[0]==TRANSVERSE_ISOTROPIC) {

            Real youngModulusTransverse = youngModulus;
            Real poissonRatioTransverse = poissonRatio;
            Real youngModulusLongitudinal = anisotropyParameter[1];

//            Real poissonRatioTransverseLongitudinal = 0.4161;
            Real poissonRatioTransverseLongitudinal = poissonRatioTransverse * sqrt(youngModulusTransverse/youngModulusLongitudinal);
            Real poissonRatioLongitudinalTransverse= poissonRatioTransverseLongitudinal*youngModulusLongitudinal/youngModulusTransverse;

            Real gamma=1/((1+poissonRatioTransverse)*(1+poissonRatioTransverse)*(1-2*poissonRatioTransverse));

            Real shearModulusLongitudinal = sqrt(youngModulusTransverse*youngModulusLongitudinal) * 1/(2*(1+poissonRatioTransverse));


            Real c11=youngModulusTransverse*gamma*(1-poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal); // = c22
            Real c33=youngModulusLongitudinal*gamma*(1-poissonRatioTransverse*poissonRatioTransverse);
            Real c44= youngModulusTransverse/(2*(1+poissonRatioTransverse)); //shearModulusLongitudinal; //
            Real c55= c44; // shearModulusLongitudinal; // = c66
            Real c12=youngModulusTransverse*gamma*(poissonRatioTransverse+poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);
            Real c13=youngModulusTransverse*gamma*(poissonRatioLongitudinalTransverse+poissonRatioTransverse*poissonRatioTransverseLongitudinal);

            // (4.119)
            Real talpha=sqrt(2.0f)*(c11+c12-c33)/(4*c13);
            Real alpha=atan(talpha);
            Real salpha=sin(alpha);
            Real calpha=cos(alpha);
            Real secalpha=1/calpha;

            // (4.118)
            Real eigen1 = c33 + M_SQRT2 * c13 * (talpha + secalpha);
            Real eigen2 = c11 - c12;
            Real eigen3 = c33 + M_SQRT2 * c13 * (talpha - secalpha);
            Real eigen4 = 2  * shearModulusLongitudinal;
///////////////////////////////////////////////////////////////////////////////////////////////////


            // ----------------------------------------------------------------------
            // BUILD THE ORTHOGONAL MATRICES
            // ----------------------------------------------------------------------

            // No dilatation/extension Kelvin modes ? Nd/Ne
            // yes see p48 :
            // because we represent the tetragonal tensor in the e3-direction
            // The five projection tensors are expressed as functions of the typical Kelvin modes and two dilatation

            // (4.120)
            //P1 = Nh1 x Nh1
            //P2 = Np3 x Np3 + Ns3 x Ns3
            //P3 = Nh2 x Nh2
            //P4 = Ns1 x Ns1 + Ns2 x Ns2

            // ---------------------------------------------------------
            Mat3x3 Nh1,Nh2;
            if (abs(poissonRatioTransverse - poissonRatioTransverseLongitudinal)> 1e-10 or abs(youngModulusLongitudinal - youngModulusTransverse) > 1e-10)
            {
                // dilatation modes
                // defined in Equation 4.121 & 4.122 (p51)
                Real val1_Nh1 = 0.5*(1+salpha)+sqrt(2.0)*calpha/4.0f;
                Real val2_Nh1 = 0.5*(1-salpha)+sqrt(2.0)*calpha/2.0f;
                Real val1_Nh2 = 0.5*(1-salpha)-sqrt(2.0)*calpha/4.0f;
                Real val2_Nh2 = 0.5*(1+salpha)-sqrt(2.0)*calpha/2.0f;

                Nh1 = val1_Nh1 * ( defaulttype::dyad(v1,v1) + defaulttype::dyad(v2,v2) ) + val2_Nh1 * defaulttype::dyad(n,n);
                Nh2 = val1_Nh2 * ( defaulttype::dyad(v1,v1) + defaulttype::dyad(v2,v2) ) + val2_Nh2 * defaulttype::dyad(n,n);

                //normalization Nh1/Nh2
                Nh1/=sqrt(2*val1_Nh1*val1_Nh1+val2_Nh1*val2_Nh1);
                Nh2/=sqrt(2*val1_Nh2*val1_Nh2+val2_Nh2*val2_Nh2);
            }
            else
            {
                if (eltIndex == 0) msg_info() << "---------> ISOTROPIC";

                Nh1; // == Nd
                Nh1.identity();
                Nh1/= sqrt(3.0f);

                Nh2=(2*defaulttype::dyad(n,n)-defaulttype::dyad(v1,v1)-defaulttype::dyad(v2,v2))/sqrt(6.0f); // == Ne3
            }

            // ---------------------------------------------------------
            // isochoric pure shear modes along e3
            // defined in Equation 4.18 (p35)
            Mat3x3 Np=(defaulttype::dyad(v1,v1)-defaulttype::dyad(v2,v2))/sqrt(2.0f);

            // ---------------------------------------------------------
            // The three isochogammaric simple shear modes along e1, e2, and e3
            // defined in Equation 4.19, Equation 4.20, and Equation 4.21, respectively (p35)
            Mat3x3 Ns1=(defaulttype::dyad(v2,n)+defaulttype::dyad(n,v2))/sqrt(2.0f);
            Mat3x3 Ns2=(defaulttype::dyad(v1,n)+defaulttype::dyad(n,v1))/sqrt(2.0f);
            Mat3x3 Ns3=(defaulttype::dyad(v1,v2)+defaulttype::dyad(v2,v1))/sqrt(2.0f);


            std::vector<Mat3x3> tmp_anisotropyMatrixArray;
            std::vector<Real> tmp_anisotropyScalarArray;


            // push all symmetric matrices and the eigenvalues
            tmp_anisotropyMatrixArray.push_back(Nh1);
            tmp_anisotropyScalarArray.push_back(eigen1);
            tmp_anisotropyMatrixArray.push_back(Np);
            tmp_anisotropyScalarArray.push_back(eigen2);
            tmp_anisotropyMatrixArray.push_back(Ns3);
            tmp_anisotropyScalarArray.push_back(eigen2);
            tmp_anisotropyMatrixArray.push_back(Nh2);
            tmp_anisotropyScalarArray.push_back(eigen3);
            tmp_anisotropyMatrixArray.push_back(Ns1);
            tmp_anisotropyScalarArray.push_back(eigen4);
            tmp_anisotropyMatrixArray.push_back(Ns2);
            tmp_anisotropyScalarArray.push_back(eigen4);


            vecAnisotropyMatrixArray[eltIndex] = tmp_anisotropyMatrixArray;
            vecAnisotropyScalarArray[eltIndex] = tmp_anisotropyScalarArray;

//            if (eltIndex == 0)
//            {
//                msg_info() << v1 << ' ' << v2 << ' '<< n;
//                msg_info() << d_anisotropyDirection.getValue()[0];
////                msg_info() << eigen1 << ' ' << eigen2 << ' '<< eigen2 << ' '<< eigen3 << ' '<< eigen4 << ' '<< eigen4 << ' '<< talpha << ' ' << secalpha;
////                msg_info() << "Nh1 " << Nh1;
////                msg_info() << "Nh2 " << Nh2;
////                msg_info() << "Np " << Np;
////                msg_info() << "Ns1 " << Ns1;
////                msg_info() << "Ns2 " << Ns2;
////                msg_info() << "Ns3 " << Ns3;
//                msg_info() << "tmp_anisotropyMatrixArray  : " << tmp_anisotropyMatrixArray[0];
//                msg_info() << "tmp_anisotropyScalarArray  : " << tmp_anisotropyScalarArray[0];
//            }
        }

    }
//    msg_info() << "-------------------------------- EXIT computeKelvinModesForElts";

}


template<class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::computeTetrahedronStiffnessEdgeMatrixForElts(size_t eltIndex,const Coord point[4], Mat6x9 edgeStiffnessVectorized[2])
{
    helper::system::thread::ctime_t startComputeStiffness=helper::system::thread::CTime::getTime();

    Coord shapeVector[4];
    Mat3x3 edgeStiffness[6];
    /// compute 6 times the rest volume
    Real volume=dot(cross(point[1]-point[0],point[2]-point[0]),point[0]-point[3]);
    /// store the rest volume
    // my_tinfo.restVolume=volume/6;

    size_t j,k,l,m,n;
    // store shape vectors at the rest configuration

    for(j=0; j<4; ++j)
    {
        if ((j%2)==0)
            shapeVector[j]=cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
        else
            shapeVector[j]= -cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;

    }

    if (elasticitySymmetry==ISOTROPIC) {
        Real mu=getMu()*fabs(volume)/6;
        Real lambda=getLambda()*fabs(volume)/6;
        Real val;

        /// compute the edge stiffness of the linear elastic material
        for(j=0; j<6; ++j)
        {
            k=edgesInTetrahedronArray[j][0];
            l=edgesInTetrahedronArray[j][1];
            // the linear stiffness matrix using shape vectors and Lame coefficients
            val=mu*dot(shapeVector[l],shapeVector[k]);
            for(m=0; m<3; ++m)
            {
                for(n=0; n<3; ++n)
                {
                    edgeStiffness[j][m][n]=lambda*shapeVector[k][n]*shapeVector[l][m]+
                        mu*shapeVector[l][n]*shapeVector[k][m];

                    if (m==n)
                    {
                        edgeStiffness[j][m][m]+=(Real)val;
                    }
                }
            }
            //msg_info() << edgeStiffness[j];
        }
        //msg_info() << "------------------------------- " << eltIndex;
    }
    else {
        size_t i;
        for(j=0; j<6; ++j)
        {
            k=edgesInTetrahedronArray[j][0];
            l=edgesInTetrahedronArray[j][1];
            // the linear stiffness matrix using shape vectors and Lame coefficients
            Mat3x3 tmp=defaulttype::dyad(shapeVector[l],shapeVector[k]);
            for(i=0;i<vecAnisotropyScalarArray[eltIndex].size();++i) {
                edgeStiffness[j]+=(vecAnisotropyScalarArray[eltIndex][i]*vecAnisotropyMatrixArray[eltIndex][i]*tmp*vecAnisotropyMatrixArray[eltIndex][i])*fabs(volume)/6;
            }
        }
    }

    size_t p;
    for(j=0; j<6; ++j)
    {
        k=edgesInTetrahedronArray[j][0];
        l=edgesInTetrahedronArray[j][1];
        for(p=0,m=0; m<3; ++m)
        {
            for(n=0; n<3; ++n,++p)
            {
                edgeStiffnessVectorized[0][j][p]=edgeStiffness[j][m][n];
                edgeStiffnessVectorized[1][j][p]=edgeStiffness[j][n][m];
            }
        }
    }
    if (this->f_printLog.getValue()) {
        helper::system::thread::ctime_t endComputeStiffness=helper::system::thread::CTime::getTime();
        totalComputeLocalStiffness+=endComputeStiffness-startComputeStiffness;
    }
}


template<class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::computeQRRotation( Mat3x3 &r, const Coord *dp)
{
    // first vector on first edge
    // second vector in the plane of the two first edges
    // third vector orthogonal to first and second

    Coord edgex = dp[0];
    edgex.normalize();

    Coord edgey = dp[1];

    Coord edgez = cross( edgex, edgey );
    edgez.normalize();

    edgey = cross( edgez, edgex );
    edgey.normalize();

    r[0][0] = edgex[0];
    r[0][1] = edgex[1];
    r[0][2] = edgex[2];
    r[1][0] = edgey[0];
    r[1][1] = edgey[1];
    r[1][2] = edgey[2];
    r[2][0] = edgez[0];
    r[2][1] = edgez[1];
    r[2][2] = edgez[2];
}


template <class DataTypes>
const helper::vector<typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::Mat3x3> &
    HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::getStiffnessArray(
    const size_t i,
    typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{
    return(restTetra->stiffnessVector);
}

template <class DataTypes>
const  helper::vector<typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::Mat3x3> &
    HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::getRotatedStiffnessArray(
    const size_t i,
    const typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{

    if (decompositionMethod==LINEAR_ELASTIC)
    {
        return(restTetra->stiffnessVector);
    } else
        return(restTetra->rotatedStiffnessVector);
}

template <class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, 
                                                                     DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & /*dataV*/ )
{
//    msg_info() << "ENTER ADDFORCE ----------------------------";
//    msg_info() << "vecAnisotropyMatrixArray  : " << vecAnisotropyMatrixArray[0][0];
//    msg_info() << "vecAnisotropyScalarArray  : " << vecAnisotropyScalarArray[0][0];

    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& x  =   dataX.getValue()  ;
    const VecCoord& x0= this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
    sofa::helper::vector<Coord> dp,force;
    size_t i,j,k,l,v0,v1,rank,p;
    size_t nbTetrahedra=_topology->getNbTetrahedra();
    
    HighOrderDegreeType degree=highOrderTetraGeo->getTopologyContainer()->getDegree();
    size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
    HighOrderTetrahedronSetTopologyContainer::VecPointID indexArray;

    if (updateTopologyInfo)
    {
        updateTopologyInformation();
    }
    helper::vector<TetrahedronRestInformation>& tetrahedronIn   f = *(tetrahedronInfo.beginEdit());
    TetrahedronRestInformation *tetinfo;
    
    dp.resize(nbControlPoints);
    force.resize(nbControlPoints);
    
    Coord dpos,sv;
        
    for(i=0; i<nbTetrahedra; i++ )
    {
        tetinfo=&tetrahedronInf[i];
        const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
            highOrderTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
        
        nbControlPoints=indexArray.size();

        if (d_oneRotationPerIntegrationPoint.getValue()) {
            
            Mat3x3 S,R;
            size_t j,k,l,m,n;
//          helper::vector<Mat3x3> stiffnessArray(nbControlPoints*(nbControlPoints-1)/2);
            helper::vector<Mat3x3> &stiffnessArray=tetinfo->rotatedStiffnessVector;
            std::fill(stiffnessArray.begin(),stiffnessArray.end(),Mat3x3());
            assert(stiffnessArray.size()==nbControlPoints*(nbControlPoints-1)/2);
            // loop through the integration points
            for (l=0;l<numericalIntegrationStiffnessDataArray.size();++l) {
                Coord dpp[6],point[4];
                // the barycentric coordinate
                highOrderTetraGeo->computeNodalValueDerivatives(i,numericalIntegrationStiffnessDataArray[l].integrationPoint, x,point);

                for(j=0; j<6; ++j){
                    m=edgesInTetrahedronArray[j][0];
                    n=edgesInTetrahedronArray[j][1];

                    dpp[j]=point[n]-point[m];

                }
                if (decompositionMethod==QR_DECOMPOSITION)
                {
                    /// perform QR decomposition
                    computeQRRotation(S,dpp);
                    R=S.transposed()*tetinfo->integrationPointsRestRotationArray[l];

                } else if (decompositionMethod==POLAR_DECOMPOSITION_MODIFIED) {

                    S[0]=dpp[0];
                    S[1]=dpp[1];
                    S[2]=dpp[2];
                    helper::Decompose<Real>::polarDecomposition( S, R );
                    R=R.transposed()*tetinfo->integrationPointsRestRotationArray[l];
                } else {
                    R.identity();
                }
                const std::vector<Vec6> & weightArray=numericalIntegrationStiffnessDataArray[l].weightArray;
                const Mat6x9  &edgeStiff1 =tetinfo->integrationPointsStiffnessVector[2*l];
                const Mat6x9  &edgeStiff2 =tetinfo->integrationPointsStiffnessVector[2*l+1];

                for (rank=0,j=0; j<nbControlPoints; ++j) {
                    v0 = indexArray[j];
                    for ( k=j+1; k<nbControlPoints; ++k,++rank) {

                        const Vec6  & coeffVec1=weightArray[2*rank];
                        const Vec6  & coeffVec2=weightArray[2*rank+1];
                        Vec9 res=edgeStiff1.multTranspose(coeffVec1)+edgeStiff2.multTranspose(coeffVec2);
                        Mat3x3 stiffness= Mat3x3((const Real *) &res[0]);
                        //          Vec9 res=edgeStiffness[0]*coeffVec1+edgeStiffness[1]*coeffVec2;


                        // loop through the integration points
                        v1 = indexArray[k];
                        dpos=x[v0]-x[v1];
                        // displacement in the rest configuration
                        dpos=R.transposed()*dpos-(x0[v0]-x0[v1]);
                        // force on first vertex in the rest configuration
                        force[k]-=R*stiffness*dpos;
                        // force on second vertex in the rest configuration
                        force[j]+=R*stiffness.multTranspose(dpos);
                        if (mparams->implicit()) {
                            // implicit scheme : need to store the rotated tensor
                            Mat3x3 mat=R*stiffness*R.transposed();
                            stiffnessArray[rank]+=mat;
                        }

                    }

                }
            }
//  tetinfo->rotatedStiffnessVector.clear();
//  tetinfo->rotatedStiffnessVector=stiffnessArray;
            for (j=0; j<nbControlPoints; ++j) {
                f[indexArray[j]]+=R*force[j];
            }


        } else {
            const  helper::vector<Mat3x3> &stiffnessArray=getStiffnessArray(i,tetinfo);
            assert(stiffnessArray.size()==nbControlPoints*(nbControlPoints-1)/2);

            if (decompositionMethod==LINEAR_ELASTIC)
            {

                for (j=0; j<nbControlPoints; ++j)
                {
                    dp[j]=x[indexArray[j]]-x0[indexArray[j]];
                }
                // loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints-1)/2
                for (rank=0,j=0; j<nbControlPoints; ++j) {
                    v0 = indexArray[j];
                    for ( k=j+1; k<nbControlPoints; ++k,++rank) {
                        v1 = indexArray[k];
                        dpos=dp[j]-dp[k];
                        //      if ((i==0))
                        //              std::cerr<<"stiffness["<<j<<","<<k<<"]="<<tetinfo->stiffnessVector[rank]<<" v0="<<v0<<" v1="<<v1<<std::endl;
                        //                   f[v1]-=tetinfo->stiffnessVector[rank]*dp[j];
                        //                   f[v0]-=tetinfo->stiffnessVector[rank].multTranspose(dp[k]);
                        f[v1]-=stiffnessArray[rank]*dpos;
                        f[v0]+=stiffnessArray[rank].multTranspose(dpos);
                    }
                }
            }
            else
            {
                Mat3x3 deformationGradient,S,R;
                Coord dpp[6];
                for (j=0; j<6; ++j)
                {
                    dpp[j]=x[tetinfo->v[edgesInTetrahedronArray[j][1]]]-x[tetinfo->v[edgesInTetrahedronArray[j][0]]];
                }
                if (decompositionMethod==POLAR_DECOMPOSITION)
                {
                    // compute the deformation gradient
                    // deformation gradient = sum of tensor product between vertex position and shape vector
                    // optimize by using displacement with first vertex
                    sv=tetinfo->shapeVector[1];


                    for (k=0; k<3; ++k)
                    {
                        for (l=0; l<3; ++l)
                        {
                            deformationGradient[k][l]=dpp[0][k]*sv[l];
                        }
                    }
                    for (j=1; j<3; ++j)
                    {
                        sv=tetinfo->shapeVector[j+1];
                        for (k=0; k<3; ++k)
                        {
                            for (l=0; l<3; ++l)
                            {
                                deformationGradient[k][l]+=dpp[j][k]*sv[l];
                            }
                        }
                    }
                    // polar decomposition of the transformation
                    helper::Decompose<Real>::polarDecomposition(deformationGradient,R);
                }
                else if (decompositionMethod==QR_DECOMPOSITION)
                {

                    /// perform QR decomposition
                    computeQRRotation(S,dpp);
                    R=S.transposed()*tetinfo->restRotation;

                } else if (decompositionMethod==POLAR_DECOMPOSITION_MODIFIED) {

                    S[0]=dpp[0];
                    S[1]=dpp[1];
                    S[2]=dpp[2];
                    helper::Decompose<Real>::polarDecomposition( S, R );
                    R=R.transposed()*tetinfo->restRotation;
                } 

                //  std::cerr<<"rotation= "<<R<<std::endl;
                //  R.identity();
                // store transpose of rotation
                tetinfo->rotation=R.transposed();
                std::fill(force.begin(),force.end(),Coord());
                if (mparams){
                    if (mparams->implicit()) {
                        tetinfo->rotatedStiffnessVector.clear();
                    }
                }
                // loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints-1)/2
                for (rank=0,j=0; j<nbControlPoints; ++j) {
                    v0 = indexArray[j];
                    for ( k=j+1; k<nbControlPoints; ++k,++rank) {
                        v1 = indexArray[k];
                        dpos=x[v0]-x[v1];
                        // displacement in the rest configuration
                        dpos=tetinfo->rotation*dpos-(x0[v0]-x0[v1]);

                        // force on first vertex in the rest configuration
                        force[k]-=stiffnessArray[rank]*dpos;
                        // force on second vertex in the rest configuration
                        force[j]+=stiffnessArray[rank].multTranspose(dpos);
                        if (mparams){
                            if (mparams->implicit()) {
                                // implicit scheme : need to store the rotated tensor
                                Mat3x3 mat=R*stiffnessArray[rank]*tetinfo->rotation;
                                tetinfo->rotatedStiffnessVector.push_back(mat);
                            }
                        }
                    }
                }
                for (j=0; j<nbControlPoints; ++j) {
                    f[indexArray[j]]+=R*force[j];
                }

            }
        }

    }

    updateMatrix=true; // next time assemble the matrix
    tetrahedronInfo.endEdit();

    dataF.endEdit();
    const double long test = dataF.getValue()[0][0];
    const double long test2 = dataF.getValue()[33][0]   ;

//    msg_info() << test <<" "<< test2;
//    msg_info() << dataF.getValue();
//    msg_info() << "EXIT ADDFORCE ----------------------------";
}


template <class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX )
{
    VecDeriv& df       = *(datadF.beginEdit());
    const VecCoord& dx =   datadX.getValue()  ;
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    const helper::vector<TetrahedronRestInformation>& tetrahedronInf = tetrahedronInfo.getValue();
    Coord dpos;
    size_t i,j,k,v0,v1,rank;
    size_t nbTetrahedra=_topology->getNbTetrahedra();
    
    HighOrderDegreeType degree=highOrderTetraGeo->getTopologyContainer()->getDegree();
    size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
    
    const TetrahedronRestInformation *tetinfo;
    

    for(i=0; i<nbTetrahedra; i++ )
    {
        tetinfo=&tetrahedronInf[i];
        const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
            highOrderTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
        const  helper::vector<Mat3x3> &stiffnessArray=getRotatedStiffnessArray(i,tetinfo);
        // create a local buffer to limit access to the df array
        sofa::helper::vector<Deriv> dforce;

        nbControlPoints=indexArray.size();
        dforce.resize(nbControlPoints);
        assert(stiffnessArray.size()==nbControlPoints*(nbControlPoints-1)/2);
        // loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
        for (rank=0,j=0; j<nbControlPoints; ++j) {
            v0 = indexArray[j];
            for ( k=j+1; k<nbControlPoints; ++k,++rank) {
                v1 = indexArray[k];
                dpos=dx[v0]-dx[v1];
                dforce[k]-=stiffnessArray[rank]*dpos*kFactor;
                dforce[j]+=stiffnessArray[rank].multTranspose(dpos*kFactor);
            }
        }
        for (j=0; j<nbControlPoints; ++j) {
            df[indexArray[j]]+=dforce[j];
        }
        /*
        indexArray.clear();
        /// get the global index of each control point in the tetrahedron
        highOrderTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfHighOrderPointsInTetrahedron(i,indexArray) ;
        if (decompositionMethod==LINEAR_ELASTIC)
        {
            // loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
            for (rank=0,j=0; j<nbControlPoints; ++j) {
                v0 = indexArray[j];
                for ( k=j+1; k<nbControlPoints; ++k,++rank) {
                    v1 = indexArray[k];
                    dpos=dx[v0]-dx[v1];
                
                    df[v1]-=tetinfo->stiffnessVector[rank]*dpos*kFactor;
                    df[v0]+=tetinfo->stiffnessVector[rank].multTranspose(dpos*kFactor);
                }
            }
        } else {
            Mat3x3 rot=tetinfo->rotation;
            // loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
            for (rank=0,j=0; j<nbControlPoints; ++j) {
                v0 = indexArray[j];
                for ( k=j+1; k<nbControlPoints; ++k,++rank) {
                    v1 = indexArray[k];
                    dpos=dx[v0]-dx[v1];
                    df[v1]-=rot.multTranspose(tetinfo->stiffnessVector[rank]*rot*dpos*kFactor);
                    df[v0]+=rot.multTranspose(tetinfo->stiffnessVector[rank].multTranspose(rot*dpos*kFactor));
                }
            }

        }
        */
    }

    datadF.endEdit();
}

template<class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix )
{
    sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    if (r)
        addKToMatrix(r.matrix, mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue()), r.offset);
    else dmsg_error() << "The function addKToMatrix found no valid matrix accessor." ;
}

template<class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal k, unsigned int &offset)
{

    const helper::vector<TetrahedronRestInformation>& tetrahedronInf = tetrahedronInfo.getValue();
    size_t i,j,l,rank,n,m;
    size_t nbTetrahedra=_topology->getNbTetrahedra();

    HighOrderDegreeType degree=highOrderTetraGeo->getTopologyContainer()->getDegree();
    size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;

    const TetrahedronRestInformation *tetinfo;

    for(i=0; i<nbTetrahedra; i++ )
    {
        tetinfo=&tetrahedronInf[i];
        const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
            highOrderTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
        const  helper::vector<Mat3x3> &stiffnessArray=getRotatedStiffnessArray(i,tetinfo);

        nbControlPoints=indexArray.size();
        assert(stiffnessArray.size()==nbControlPoints*(nbControlPoints-1)/2);

        // loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
        for (rank=0,j=0; j<nbControlPoints; ++j) {
            for ( l=j+1; l<nbControlPoints; ++l,++rank) {
                for (n=0; n<3; ++n) {
                    for (m=0; m<3; ++m) {
                        Mat3x3 stiffnessArrayTranspose = stiffnessArray[rank]   ;
                        stiffnessArrayTranspose.transpose();
                        mat->add(3*indexArray[l]+ n + offset, 3*indexArray[l]+ m + offset, stiffnessArray[rank][n][m] * k);
                        mat->add(3*indexArray[l]+ n + offset, 3*indexArray[j]+ m + offset, - stiffnessArray[rank][n][m] * k);
                        mat->add(3*indexArray[j]+ n + offset, 3*indexArray[l]+ m + offset, - stiffnessArrayTranspose[n][m] * k);
                        mat->add(3*indexArray[j]+ n + offset, 3*indexArray[j]+ m + offset, stiffnessArrayTranspose[n][m] * k);
                    }
                }
            }
        }
    }
}
