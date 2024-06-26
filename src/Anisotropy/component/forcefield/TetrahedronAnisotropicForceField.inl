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

#include <Anisotropy/component/forcefield/TetrahedronAnisotropicForceField.h>
#include <sofa/core/topology/TopologyData.inl>

#include <sofa/helper/ColorMap.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/gl/gl.h>


namespace anisotropy::forcefield
{

const unsigned int edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};


template <class DataTypes>
TetrahedronAnisotropicForceField<DataTypes>::TetrahedronAnisotropicForceField()
    : tetrahedronInfo(initData(&tetrahedronInfo, "tetrahedronInfo", "Internal tetrahedron data"))
    , updateMatrix(true)
    , d_poissonRatio(initData(&d_poissonRatio, "poissonRatio", "Poisson ratio in Hooke's law"))
    , d_youngModulus(initData(&d_youngModulus, "youngModulus", "Young modulus in Hooke's law"))
    , d_anisotropy(initData(&d_anisotropy, std::string("isotropic"), "elasticitySymmetry", "the type of anisotropy for the elasticity tensor :\"isotropic\"  or \"transverseIsotropic\" or \"ortho_scalar\" or \"cubic\" "))
    , d_anisotropyParameter(initData(&d_anisotropyParameter, "anisotropyParameters", "the elastic parameters for anisotropic materials.\n"
                                                                                     "- for cubic symmetry       --> anisotropyParameters == [anisotropyRatio]\n"
                                                                                     "- for transverse symemetry --> anisotropyParameters == [youngModulusLongitudinal, poissonRatioTransverseLongitudinal, shearModulusTransverse]"))
    , d_anisotropyDirection(initData(&d_anisotropyDirection, "anisotropyDirections", "the directions of anisotropy"))
    , d_controlPoints(initData(&d_controlPoints,"controlPoints","controlPoints"))
    , d_IDWDepth(initData(&d_IDWDepth,2,"IDWDepth","How many CP a data is interpolated upon"))
    , d_meshRotation(initData(&d_meshRotation,"meshRotation",""))
    , d_drawHeterogeneousTetra(initData(&d_drawHeterogeneousTetra,true,"drawHeterogeneousTetra","Draw Heterogeneous Tetra in different color"))
    , d_drawDirection(initData(&d_drawDirection,true,"drawDirection","Draw different color for each direction"))
    , d_transparency(initData(&d_transparency,(Real)0.25,"transparency","transparency [0,1]"))
    , lambda(0)
    , mu(0)
{
//numericalIntegrationMethod='Tetrahedron Gauss'
//integrationMethod='analytical'
//method="qr"
//integrationOrder=1
//orderDegree=1,
}


template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);

    this->Inherited::init();


    if (d_anisotropy.getValue() == "transverseIsotropic")
        elasticitySymmetry= TRANSVERSE_ISOTROPIC;
    else
        elasticitySymmetry= ISOTROPIC;

    /// Init of tetrahedronInf which is a container which has, among other, mechanical info for each tetra:
    /// shapeVector, rotation, stiffnessVector ...
    _topology = this->getContext()->getMeshTopology();
    const std::vector< Tetrahedron > &tetrahedronArray = this->_topology->getTetrahedra();
    type::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    tetrahedronInf.resize(_topology->getNbTetrahedra());

    /// We are putting an elasticity tensor (E) for each tetra of our topology
    /// This 4d order elasticity tensor will be projected thanks to Kelvin modes
    /// ref : On Inverse Form Finding for Anisotropic Materials in the Logarithmic Strain Space, Sandrine Germain, 2013
    if (d_anisotropyDirection.getValue().size() == _topology->getNbTetrahedra() or d_controlPoints.isSet())
    {
        /// Each tetra E will be decompose thanks to its modes.
        /// This decomposition will result in :
        /// - 6 eigens values (lambda)
        /// - 6 eigentensors (3x3 matrix) which when using the Kronecker product on themselves will give the projection operators (P).
        /// Then we have this decomposition: E = sum([k=1->N_modes] of lambda_k * P_k) (4.3)
        for (size_t i=0; i<_topology->getNbTetrahedra(); ++i)
        {
            if (d_controlPoints.getValue().size() > 0)
            {
                setMechanicalParametersFromControlPoints(i,tetrahedronArray[i]);
                // msg_info() << d_anisotropyParameter.getValue()[i];
                // msg_info() << d_anisotropyDirection.getValue()[i];
            }
            /// Given as entries, mechanical values (up to 9 if Orthotropic):
            /// (Young_Modulus_(x/y/z), Poisson_ratio_(x/y/z), Shear_Modulus_(x/y/z)
            /// We initialize each tetrahedronInf with the eigen tensors/values associated
            computeKelvinModesForElts(i,tetrahedronInf[i].eigenTensors,tetrahedronInf[i].eigenValues);

            /// Init tetrahedronInf restRotation & stiffnessVector
            initStiffnessVector(tetrahedronInf[i],tetrahedronArray[i]);
        }
    }
    else
        msg_error() << "anisotropyDirection dim != NbTetrahedra in mesh " << d_controlPoints.isSet();


    tetrahedronInfo.endEdit();
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    updateTopologyInfo=true;
}


template <class DataTypes>
inline void TetrahedronAnisotropicForceField<DataTypes>::reinit()
{
//    updateStiffnessVectorWithCP();
}

template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::updateTopologyInformation()
{
    int i;
    unsigned int j;

    int nbTetrahedra=_topology->getNbTetrahedra();

    TetrahedronRestInformation *tetinfo;

    type::vector<typename TetrahedronAnisotropicForceField<DataTypes>::TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());

    for(i=0; i<nbTetrahedra; i++ )
    {
        tetinfo=&tetrahedronInf[i];
        /// describe the jth vertex index of triangle no i
        const Tetrahedron &ta= _topology->getTetrahedron(i);

        for (j=0; j<4; ++j)
        {
            tetinfo->v[j]=ta[j];
        }
    }
    updateTopologyInfo=false;
    tetrahedronInfo.endEdit();
}

template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::setMechanicalParametersFromControlPoints(size_t eltIndex,Tetra indexArray)
{
    // msg_info() << "IDW ELT " << eltIndex;

    AnisotropyDirectionArray &dir = *d_anisotropyDirection.beginEdit();
    type::vector<ParameterArray> &param = *d_anisotropyParameter.beginEdit();
    type::vector<Real> &youngModulus = *d_youngModulus.beginEdit();
    type::vector<Real> &poissonRatio = *d_poissonRatio.beginEdit();

    if (updateTopologyInfo and eltIndex == 0)
    {
        dir.resize(_topology->getNbTetrahedra());
        param.resize(_topology->getNbTetrahedra());
        youngModulus.resize(_topology->getNbTetrahedra());
        poissonRatio.resize(_topology->getNbTetrahedra());
        for (size_t i=0; i<_topology->getNbTetrahedra(); ++i)
        {
            param[i] = {2,0.,0.4161,0.};
            poissonRatio[i] = 0.4161; // TEMPORARY HARD SET
            // param[i] = {2,0.,0.26162,0.};
            // poissonRatio[i] = 0.26162; // TEMPORARY HARD SET
        }
    }

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
    const Coord tetraBarycenter= (x[indexArray[0]] + x[indexArray[1]] + x[indexArray[2]] + x[indexArray[3]]) / 4;

    type::vector<std::pair<Real, int> > listOfNorm = generateListOfNorm(tetraBarycenter);

    IDWdata(listOfNorm,youngModulus[eltIndex],1); // Young Modulus Transverse
    IDWdata(listOfNorm,param[eltIndex][1],2);     // Young Modulus Longituinal
    IDWdirection(listOfNorm,dir[eltIndex]);       // direction of transverse anisotropy vec normalized Coord(x,y,z)

    d_anisotropyDirection.endEdit();
    d_anisotropyParameter.endEdit();
    d_youngModulus.endEdit();
    d_poissonRatio.endEdit();
}

template <class DataTypes>
type::vector<std::pair< typename TetrahedronAnisotropicForceField<DataTypes>::Real, int>> TetrahedronAnisotropicForceField<DataTypes>::generateListOfNorm(const Coord tetraBarycenter)
{
    const type::vector<Vec5>&  controlPoints = d_controlPoints.getValue();
    type::vector<Real> listOfNorm;

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();

    for (size_t j=0;j<controlPoints.size();j++)
    {
        listOfNorm.push_back((tetraBarycenter-x[controlPoints[j][0]]).norm());
    }

    type::vector<std::pair<Real, int> > vp;
    sortArr(listOfNorm, listOfNorm.size(),&vp);

    Real sumWeights = 0;
    size_t power = 3;
    for(Index i=0;i<d_IDWDepth.getValue();i++)
    {
        vp[i].first = 1/(pow(vp[i].first,power));
        sumWeights += vp[i].first;
    }

    for (Index i=0;i<d_IDWDepth.getValue();i++)
    {
        vp[i].first /= sumWeights;
    }

    return vp;
}

template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::IDWdata(type::vector<std::pair<Real, int> > listOfNormWeighted, Real& dataToInterpolate, size_t indexData)
{
    Real &data = dataToInterpolate;
    const type::vector<Vec5>&  controlPoints = d_controlPoints.getValue();
    // msg_info() << "indexData " << indexData << "    " << controlPoints[listOfNormWeighted[0].second][indexData];
    for(size_t i=0;i<d_IDWDepth.getValue();i++)
    {
        data += controlPoints[listOfNormWeighted[i].second][indexData] * listOfNormWeighted[i].first;
    }
}

template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::IDWdirection(type::vector<std::pair<Real, int> > listOfNormWeighted, Coord &dir)
{
    const type::vector<Vec5>&  controlPoints = d_controlPoints.getValue();
    const Real meshRotation = d_meshRotation.getValue();

    type::vector<Coord> listOfDir;
    for (size_t j=0;j<controlPoints.size();j++)
    {
        /// DON'T FORGET MESH ROTATION FOR BEAMS
        //listOfDir.push_back(Coord(cos(controlPoints[j][3])*cos(controlPoints[j][4]+meshRotation),sin(controlPoints[j][3])*cos(controlPoints[j][4]+meshRotation),sin(controlPoints[j][4]+meshRotation)));
        /// ELSE
        //listOfDir.push_back(Coord(cos(controlPoints[j][3])*cos(controlPoints[j][4]),sin(controlPoints[j][3])*cos(controlPoints[j][4]),sin(controlPoints[j][4])));
        Coord tmp = Coord(cos(controlPoints[j][3])*cos(controlPoints[j][4]),sin(controlPoints[j][3])*cos(controlPoints[j][4]),sin(controlPoints[j][4]));
        /// MESH ROTATION AROUND X AXIS
        Real y = tmp[1];
        tmp[1] = tmp[1]*cos(meshRotation)-tmp[2]*sin(meshRotation);
        tmp[2] = y*sin(meshRotation)+tmp[2]*cos(meshRotation);
        listOfDir.push_back(tmp);
    }

    Real d_x=0,d_y=0,d_z=0;

    for(size_t i=0;i<d_IDWDepth.getValue();i++)
    {
        d_x += listOfDir[listOfNormWeighted[i].second][0] * listOfNormWeighted[i].first;
        d_y += listOfDir[listOfNormWeighted[i].second][1] * listOfNormWeighted[i].first;
        d_z += listOfDir[listOfNormWeighted[i].second][2] * listOfNormWeighted[i].first;
    }
    dir = Coord(d_x,d_y,d_z);
}

template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::computeKelvinModesForElts(size_t eltIndex,EigenTensors &eigenTensors , EigenValues &eigenValues)
{
    // msg_info() << "-------------------------------- ENTER computeKelvinModesForElts";

    if (elasticitySymmetry != ISOTROPIC)
    {
        Vec4 anisotropyParameter=d_anisotropyParameter.getValue()[eltIndex];
        // Vec10 anisotropyParameter=d_anisotropyParameter.getValue()[eltIndex];
        Real youngModulus = d_youngModulus.getValue()[eltIndex];
        Real poissonRatio = d_poissonRatio.getValue()[eltIndex];

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

        /// Here we will use the work of Sandrine Germian see thesis (2015): https://opus4.kobv.de/opus4-fau/frontdoor/index/index/docId/3490
        /// Chap 4 Spectral decomposition and the Kelvin modes :
        ///
        /// - For CUBIC                  -------> 4.4.1 The cubic crystal system Materials (p38)
        ///
        /// - For TRANSVERSE_ISOTROPIC   -------> 4.4.6 The tetragonal crystal system Materials (p48)
        ///
        /// As well as the equations from http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm

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

            Mat3x3 Ne=(2*type::dyad(n,n)-type::dyad(v1,v1)-type::dyad(v2,v2))/sqrt(6.0f);
            Mat3x3 Np=(type::dyad(v1,v1)-type::dyad(v2,v2))/sqrt(2.0f);

            // The third projection tensor associated with the third eigenvalue is the sum of the tensor product of the three simple shear modes with themselves
            // P3 = Ns1 x Ns1 + Ns2 x Ns2 + Ns3 x Ns3

            Mat3x3 Ns1=(type::dyad(v2,n)+type::dyad(n,v2))/sqrt(2.0f);
            Mat3x3 Ns2=(type::dyad(v1,n)+type::dyad(n,v1))/sqrt(2.0f);
            Mat3x3 Ns3=(type::dyad(v1,v2)+type::dyad(v2,v1))/sqrt(2.0f);


            // push all symmetric matrices and the eigenvalues
//            std::vector<Mat3x3> eigenTensors;
//            std::vector<Real> anisotropyScalarArray;

//            eigenTensors.push_back(Nd);
//            anisotropyScalarArray.push_back(eigen1);
//            eigenTensors.push_back(Ne);
//            anisotropyScalarArray.push_back(eigen2);
//            eigenTensors.push_back(Np);
//            anisotropyScalarArray.push_back(eigen2);
//            eigenTensors.push_back(Ns1);
//            anisotropyScalarArray.push_back(eigen3);
//            eigenTensors.push_back(Ns2);
//            anisotropyScalarArray.push_back(eigen3);
//            eigenTensors.push_back(Ns3);
//            anisotropyScalarArray.push_back(eigen3);

//            veceigenTensors.push_back(eigenTensors);
//            vecAnisotropyScalarArray.push_back(anisotropyScalarArray);

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
            // Real youngModulusTransverse = youngModulus;
            // Real poissonRatioTransverse = poissonRatio;
            // Real youngModulusLongitudinal = anisotropyParameter[1];

            // //            Real poissonRatioTransverseLongitudinal = 0.4161;
            // Real poissonRatioTransverseLongitudinal = poissonRatioTransverse * sqrt(youngModulusTransverse/youngModulusLongitudinal);
            // Real poissonRatioLongitudinalTransverse = poissonRatioTransverseLongitudinal*youngModulusLongitudinal/youngModulusTransverse;

            // Real gamma=1/((1+poissonRatioTransverse)*(1+poissonRatioTransverse)*(1-2*poissonRatioTransverse));

            // Real shearModulusLongitudinal = sqrt(youngModulusTransverse*youngModulusLongitudinal) * 1/(2*(1+poissonRatioTransverse));


            // Real c11=youngModulusTransverse*gamma*(1-poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal); // = c22
            // Real c33=youngModulusLongitudinal*gamma*(1-poissonRatioTransverse*poissonRatioTransverse);
            // Real c44= youngModulusTransverse/(2*(1+poissonRatioTransverse)); //shearModulusLongitudinal; //
            // Real c55= c44; // shearModulusLongitudinal; // = c66
            // Real c12=youngModulusTransverse*gamma*(poissonRatioTransverse+poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);
            // Real c13=youngModulusTransverse*gamma*(poissonRatioLongitudinalTransverse+poissonRatioTransverse*poissonRatioTransverseLongitudinal);


            ///////////////////////////////////////////////////////////////////////////////////////////////////


            // // get the constants from the young modulus, Poisson ratio and anisotropy ratio.
            // long double youngModulusTransverse = youngModulus;
            // long double poissonRatioTransverse = poissonRatio;
            // long double youngModulusLongitudinal=anisotropyParameter[1];

            // long double poissonRatioTransverseLongitudinal=anisotropyParameter[2];
            // //Real poissonRatioTransverseLongitudinal = poissonRatioTransverse * sqrt(youngModulusTransverse/youngModulusLongitudinal);

            // long double poissonRatioLongitudinalTransverse=poissonRatioTransverseLongitudinal*youngModulusLongitudinal/youngModulusTransverse;

            // //Real shearModulusLongitudinal=anisotropyParameter[3];
            // //Real shearModulusLongitudinal = youngModulusLongitudinal/(2.0*(1.0+poissonRatioTransverseLongitudinal));
            // long double shearModulusLongitudinal = youngModulusLongitudinal/(2.0*(1.0+poissonRatioLongitudinalTransverse));


            // //if (poissonRatioLongitudinalTransverse>0.5)
            // //    poissonRatioLongitudinalTransverse = 0.5;

            // long double gamma=1/(1-pow(poissonRatioTransverse,2)-2*poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal-2*poissonRatioTransverse*poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);

            // long double c11=youngModulusTransverse*gamma*(1-poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal); // = c22
            // long double c33=youngModulusLongitudinal*gamma*(1-pow(poissonRatioTransverse,2));

            // long double c44= youngModulusTransverse/(2*(1+poissonRatioTransverse)); //shearModulusLongitudinal; //
            // long double c55= shearModulusLongitudinal; //c44; // = c66
            // long double c12=youngModulusTransverse*gamma*(poissonRatioTransverse+poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);

            // // STABLE
            // long double c13=youngModulusTransverse*gamma*(poissonRatioLongitudinalTransverse+poissonRatioTransverse*poissonRatioTransverseLongitudinal);
            // // UNSTABLE BUT CORRECT FORMULATION !!?
            // // long double c13 = youngModulusLongitudinal*gamma*(poissonRatioTransverseLongitudinal+poissonRatioTransverse*poissonRatioTransverseLongitudinal);


            ///////////////////////////////////////////////////////////////////////////////////////////////////


            // get the constants from the young modulus, Poisson ratio and anisotropy ratio.
            Real youngModulusTransverse = youngModulus;
            Real poissonRatioTransverse = poissonRatio;
            Real youngModulusLongitudinal=anisotropyParameter[1];
            Real poissonRatioTransverseLongitudinal=anisotropyParameter[2];
            Real shearModulusLongitudinal=anisotropyParameter[3];
            //Real poissonRatioTransverseLongitudinal = poissonRatioTransverse * sqrt(youngModulusTransverse/youngModulusLongitudinal);
            //            Real shearModulusLongitudinal = youngModulusLongitudinal/(2.0*(1.0+poissonRatioTransverseLongitudinal));

            Real poissonRatioLongitudinalTransverse= poissonRatioTransverseLongitudinal*youngModulusLongitudinal/youngModulusTransverse;

            Real gamma=1/(1-pow(poissonRatioTransverse,2)-2*poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal-2*poissonRatioTransverse*poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);

            Real c11=youngModulusTransverse*gamma*(1-poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal); // = c22
            Real c33=youngModulusLongitudinal*gamma*(1-poissonRatioTransverse*poissonRatioTransverse);
            Real c44= youngModulusTransverse/(2*(1+poissonRatioTransverse)); //shearModulusLongitudinal; //
            Real c55= c44; // shearModulusLongitudinal; // = c66
            Real c12=youngModulusTransverse*gamma*(poissonRatioTransverse+poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);
            Real c13=youngModulusTransverse*gamma*(poissonRatioLongitudinalTransverse+poissonRatioTransverse*poissonRatioTransverseLongitudinal);


            ///////////////////////////////////////////////////////////////////////////////////////////////////



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
            Real eigen4 = 2  * c44; //shearModulusLongitudinal;
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

                Nh1 = val1_Nh1 * ( type::dyad(v1,v1) + type::dyad(v2,v2) ) + val2_Nh1 * type::dyad(n,n);
                Nh2 = val1_Nh2 * ( type::dyad(v1,v1) + type::dyad(v2,v2) ) + val2_Nh2 * type::dyad(n,n);

                //normalization Nh1/Nh2
                Nh1/=sqrt(2*val1_Nh1*val1_Nh1+val2_Nh1*val2_Nh1);
                Nh2/=sqrt(2*val1_Nh2*val1_Nh2+val2_Nh2*val2_Nh2);

            }
            else
            {
                if (eltIndex == 0) msg_info() << "---------> ISOTROPIC";
                v1 = Coord(1,0,0);
                v2 = Coord(0,1,0);
                n  = Coord(0,0,1);

                // (4.29)
                eigen1 = c11 + 2*c12;
                eigen2 = c11 - c12;
                eigen3 = eigen2;
                eigen4 = eigen2;

                Nh1; // == Nd
                Nh1.identity();
                Nh1/= sqrt(3.0f);

                Nh2=(2*type::dyad(n,n)-type::dyad(v1,v1)-type::dyad(v2,v2))/sqrt(6.0f); // == Ne3
            }

            // ---------------------------------------------------------
            // isochoric pure shear modes along e3
            // defined in Equation 4.18 (p35)
            Mat3x3 Np=(type::dyad(v1,v1)-type::dyad(v2,v2))/sqrt(2.0f);


            // ---------------------------------------------------------
            // The three isochogammaric simple shear modes along e1, e2, and e3
            // defined in Equation 4.19, Equation 4.20, and Equation 4.21, respectively (p35)
            Mat3x3 Ns1=(type::dyad(v2,n)+type::dyad(n,v2))/sqrt(2.0f);
            Mat3x3 Ns2=(type::dyad(v1,n)+type::dyad(n,v1))/sqrt(2.0f);
            Mat3x3 Ns3=(type::dyad(v1,v2)+type::dyad(v2,v1))/sqrt(2.0f);



            // push all symmetric matrices and the eigenvalues
            eigenTensors[0] = Nh1;
            eigenValues[0] = eigen1;
            eigenTensors[1] = Np;
            eigenValues[1] = eigen2;
            eigenTensors[2] = Ns3;
            eigenValues[2] = eigen2;
            eigenTensors[3] = Nh2;
            eigenValues[3] = eigen3;
            eigenTensors[4] = Ns1;
            eigenValues[4] = eigen4;
            eigenTensors[5] = Ns2;
            eigenValues[5] = eigen4;

            if (eltIndex == 0)
            {
                msg_info() << v1 << ' ' << v2 << ' '<< n;
                msg_info() << d_anisotropyDirection.getValue()[0];
                msg_info() << eigen1 << ' ' << eigen2 << ' '<< eigen2 << ' '<< eigen3 << ' '<< eigen4 << ' '<< eigen4 << ' '<< talpha << ' ' << secalpha;
                msg_info() << "Nh1 " << Nh1;
                msg_info() << "Nh2 " << Nh2;
                msg_info() << "Np " << Np;
                msg_info() << "Ns1 " << Ns1;
                msg_info() << "Ns2 " << Ns2;
                msg_info() << "Ns3 " << Ns3;
                msg_info() << "tmp_anisotropyMatrixArray  : " << eigenTensors[0];
                msg_info() << "tmp_anisotropyScalarArray  : " << eigenValues[0];
            }
        }

    }
    // msg_info() << "-------------------------------- EXIT computeKelvinModesForElts";

}

template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::initStiffnessVector(TetrahedronRestInformation &my_tinfo,
                                                                      const Tetrahedron &tt)
{
    // set array to zero
    std::fill(my_tinfo.stiffnessVector.begin(), my_tinfo.stiffnessVector.end(), Mat3x3());


    typename DataTypes::Coord point[4];
    const typename DataTypes::VecCoord &restPosition = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();

    // store the point position
    for (size_t j = 0; j < 4; ++j)
        point[j] = (restPosition)[tt[j]];

    size_t k, l;
    for (size_t  j = 0; j < 6; ++j)
    {
        k = edgesInTetrahedronArray[j][0];
        l = edgesInTetrahedronArray[j][1];

        // store the rest edge vector
        my_tinfo.restEdgeVector[j] = point[l] - point[k];
    }

    /// compute the rotation matrix of the initial tetrahedron for the QR decomposition
    /// Fill restRotation with previously computed restEdgeVector
    computeQRRotation(my_tinfo.restRotation, my_tinfo.restEdgeVector);

    /// Filling tetrahedronInf->stiffnessVector thanks to previous E decomposition
    /// Important: the stiffnessVector is here defined on the edges of tetra
    computeTetrahedronStiffnessEdgeMatrixForElts(my_tinfo.stiffnessVector,my_tinfo.eigenTensors,my_tinfo.eigenValues,point);

    // msg_info("my_tinfo.stiffnessVector") << my_tinfo.stiffnessVector;
}

template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::updateStiffnessVector(const Tetrahedron &tetra,TetrahedronRestInformation &my_tinfo)
{
    const typename DataTypes::VecCoord &position = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();

    // set array to zero
    std::fill(my_tinfo.stiffnessVector.begin(), my_tinfo.stiffnessVector.end(), Mat3x3());

    typename DataTypes::Coord point[4];
    // store the point position
    for (size_t j = 0; j < 4; ++j)
        point[j] = (position)[tetra[j]];

    computeTetrahedronStiffnessEdgeMatrixForElts(my_tinfo.stiffnessVector,my_tinfo.eigenTensors,my_tinfo.eigenValues,point);
}

template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::computeTetrahedronStiffnessEdgeMatrixForElts(TetraEdgesStiffness &stiffnessVector,const EigenTensors &eigenTensors , const EigenValues &eigenValues,const Coord point[4])
{
    Coord shapeVector[4];
    /// compute 6 times the rest volume
    Real volume=dot(cross(point[1]-point[0],point[2]-point[0]),point[0]-point[3]);

    size_t j,k,l,m,n;
    /// store shape vectors at the rest configuration
    for(j=0; j<4; ++j)
    {
        if ((j%2)==0)
            shapeVector[j]= cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
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
                    stiffnessVector[j][m][n]=lambda*shapeVector[k][n]*shapeVector[l][m]+
                                             mu*shapeVector[l][n]*shapeVector[k][m];

                    if (m==n)
                    {
                        stiffnessVector[j][m][m]+=(Real)val;
                    }
                }
            }
        }
    }
    else {
        size_t i;
        for(j=0; j<6; ++j)
        {
            k=edgesInTetrahedronArray[j][0];
            l=edgesInTetrahedronArray[j][1];

            /// Stiffness Matrix decomposition thanks the kelvin modes and shape vector
            /// ref : AnisotropicElasticity Hervé Delingette, 2019, (4)
            Mat3x3 tmp=type::dyad(shapeVector[l],shapeVector[k]);
            for(i=0;i<6;++i) {
                stiffnessVector[j]+=(eigenValues[i]*eigenTensors[i]*tmp*eigenTensors[i])*fabs(volume)/6;
            }
        }
    }
}


template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::computeQRRotation( Mat3x3 &r, const Coord *dp)
{
    /// first vector on first edge
    /// second vector in the plane of the two first edges
    /// third vector orthogonal to first and second

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
const typename  TetrahedronAnisotropicForceField<DataTypes>::TetraEdgesStiffness& TetrahedronAnisotropicForceField<DataTypes>::getStiffnessArray(
   const typename TetrahedronAnisotropicForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{
    return(restTetra->stiffnessVector);
}

template <class DataTypes>
const  typename TetrahedronAnisotropicForceField<DataTypes>::TetraEdgesStiffness&
TetrahedronAnisotropicForceField<DataTypes>::getRotatedStiffnessArray(
    const typename TetrahedronAnisotropicForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{
    return(restTetra->rotatedStiffnessVector);
}

template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */,
                                                                        DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & /*dataV*/ )
{
    // msg_info() << "ENTER addForce";

    SOFA_UNUSED(mparams);

    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& x  = dataX.getValue();
    const VecCoord& x0= this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
    sofa::type::vector<Coord> dp,force;
    size_t i,j,k,l;
    size_t v0,v1;
    const VecElement * _indexedElements = & (_topology->getTetrahedra());

    if (updateTopologyInfo)
    {
        updateTopologyInformation();
    }

    type::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    TetrahedronRestInformation *tetinfo;


    dp.resize(4);
    force.resize(4);

    Coord dpos,sv;
    typename VecElement::const_iterator it;
    for(it=_indexedElements->begin(), i = 0 ; it!=_indexedElements->end(); ++it,++i)
    {
        Tetra index = *it;
        tetinfo=&tetrahedronInf[i];
        size_t nbControlPoints=index.size();
        const  TetraEdgesStiffness &stiffnessArray=getStiffnessArray(tetinfo);


        Mat3x3 deformationGradient,S,R;
        Coord dpp[6];
        for (j=0; j<6; ++j)
        {
            dpp[j]=x[tetinfo->v[edgesInTetrahedronArray[j][1]]]-x[tetinfo->v[edgesInTetrahedronArray[j][0]]];
        }

        /// perform QR decomposition
        computeQRRotation(S,dpp);
        R=S.transposed()*tetinfo->restRotation;

        /// store transpose of rotation
        tetinfo->rotation=R.transposed();
        std::fill(force.begin(),force.end(),Coord());

        tetinfo->rotatedStiffnessVector.clear();

        // loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints-1)/2
        for (l=0,j=0; j<nbControlPoints; ++j) {
            v0 = index[j];
            for ( k=j+1; k<nbControlPoints; ++k,++l) {
                v1 = index[k];
                dpos=x[v0]-x[v1];
                // displacement in the rest configuration
                dpos=tetinfo->rotation*dpos-(x0[v0]-x0[v1]);
                // force on first vertex in the rest configuration
                force[k]-=stiffnessArray[l]*dpos;
                // force on second vertex in the rest configuration
                force[j]+=stiffnessArray[l].multTranspose(dpos);

                // implicit scheme : need to store the rotated tensor
                Mat3x3 mat=R*stiffnessArray[l]*tetinfo->rotation;
                tetinfo->rotatedStiffnessVector[l] = mat;

            }
        }

        for (j=0; j<nbControlPoints; ++j)
        {
            f[index[j]]+=R*force[j];
        }

    }
    updateMatrix=true; // useless normally
    tetrahedronInfo.endEdit();
    dataF.endEdit();
}


template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX )
{
    // msg_info() << "ENTER addDForce";

    VecDeriv& df       = *(datadF.beginEdit());
    const VecCoord& dx =   datadX.getValue()  ;
    Real kFactor = sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams,this->rayleighStiffness.getValue());
    type::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());

    Coord dpos;
    size_t i,j,k,v0,v1,rank;

    size_t nbControlPoints= 4;

    TetrahedronRestInformation *tetinfo;
    const VecElement * _indexedElements = & (_topology->getTetrahedra());


    typename VecElement::const_iterator it;
    for(it=_indexedElements->begin(), i = 0 ; it!=_indexedElements->end(); ++it,++i)
    {
        Tetra index = *it;
        tetinfo=&tetrahedronInf[i];
        const  TetraEdgesStiffness &stiffnessArray=getRotatedStiffnessArray(tetinfo);
        sofa::type::vector<Deriv> dforce;

        dforce.resize(nbControlPoints);

        /// loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
        for (rank=0,j=0; j<nbControlPoints; ++j)
        {
            v0 = index[j];
            for ( k=j+1; k<nbControlPoints; ++k,++rank)
            {
                v1 = index[k];
                dpos=dx[v0]-dx[v1];
                dforce[k]-=stiffnessArray[rank]*dpos*kFactor;
                dforce[j]+=stiffnessArray[rank].multTranspose(dpos*kFactor);
            }
        }

        for (j=0; j<nbControlPoints; ++j)
        {
            df[index[j]]+=dforce[j];
        }
    }

    datadF.endEdit();
    tetrahedronInfo.endEdit();
}


template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::addKToMatrix(sofa::linearalgebra::BaseMatrix *mat, SReal k, unsigned int &offset)
{
    msg_info() << "ENTER addKToMatrix";

    const type::vector<TetrahedronRestInformation> tetrahedronInf = *(tetrahedronInfo.beginEdit());
    size_t i,j,l,rank,n,m;

    size_t nbControlPoints=4;

    const TetrahedronRestInformation *tetinfo;
    const VecElement * _indexedElements = & (_topology->getTetrahedra());


    typename VecElement::const_iterator it;
    for(it=_indexedElements->begin(), i = 0 ; it!=_indexedElements->end(); ++it,++i)
    {
        Tetra index = *it;
        tetinfo=&tetrahedronInf[i];
        const  TetraEdgesStiffness &stiffnessArray=getStiffnessArray(tetinfo);
        // const  TetraEdgesStiffness &stiffnessArray=getRotatedStiffnessArray(tetrahedronInf[i]);

        // loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
        for (rank=0,j=0; j<nbControlPoints; ++j) {
            for ( l=j+1; l<nbControlPoints; ++l,++rank) {
                for (n=0; n<3; ++n) {
                    for (m=0; m<3; ++m) {
                        Mat3x3 stiffnessArrayTranspose = stiffnessArray[rank];
                        stiffnessArrayTranspose.transpose();
                        mat->add(3*index[l]+ n + offset, 3*index[l]+ m + offset, stiffnessArray[rank][n][m] * k);
                        mat->add(3*index[l]+ n + offset, 3*index[j]+ m + offset, - stiffnessArray[rank][n][m] * k);
                        mat->add(3*index[j]+ n + offset, 3*index[l]+ m + offset, - stiffnessArrayTranspose[n][m] * k);
                        mat->add(3*index[j]+ n + offset, 3*index[j]+ m + offset, stiffnessArrayTranspose[n][m] * k);
                    }
                }
            }
        }
    }
    tetrahedronInfo.endEdit();
}


template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::updateLameCoefficients()
{

    lambda= d_youngModulus.getValue()[0]*d_poissonRatio.getValue()[0]/((1-2*d_poissonRatio.getValue()[0])*(1+d_poissonRatio.getValue()[0])); // E*v/(1-2*v)*(1+v)
    mu = d_youngModulus.getValue()[0]/(2*(1+d_poissonRatio.getValue()[0])); // E/(2*(1+v))

}


template<class DataTypes>
SReal TetrahedronAnisotropicForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const
{
    msg_error() << "ERROR("<<this->getClassName()<<"): getPotentialEnergy( const MechanicalParams*, const DataVecCoord& ) not implemented.";
    return 0.0;
}


template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

#ifndef SOFA_NO_OPENGL
    bool test_visu = 1;
    if(test_visu == 0){
        if (!vparams->displayFlags().getShowForceFields()) return;
        if (!this->mstate) return;

        if (vparams->displayFlags().getShowWireFrame())
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);


        if (vparams->displayFlags().getShowWireFrame())
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
    else
    {
        if(this->d_componentState.getValue() == sofa::core::objectmodel::ComponentState::Invalid) return ;
        if (!vparams->displayFlags().getShowForceFields()) return;

        type::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
        const VecElement * _indexedElements = & (_topology->getTetrahedra());

        const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();

        size_t i;
        size_t nbTetrahedra=_topology->getNbTetrahedra();

        Real young_1, young_2,anisotropyType;
        float transparency = d_transparency.getValue();
        typename VecElement::const_iterator it;

        for(it=_indexedElements->begin(), i = 0 ; it!=_indexedElements->end(); ++it,++i)
        {
            anisotropyType = d_anisotropyParameter.getValue()[i][0];


            Tetra indexArray = *it;
            std::vector< type::Vec3 > points[4];

            //nbControlPoints=indexArray.size();
            //assert(stiffnessArray.size()==nbControlPoints*(nbControlPoints-1)/2);
            //// loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
            //for (rank=0,j=0; j<nbControlPoints; ++j) {

            Index a = indexArray[0];
            Index b = indexArray[1];
            Index c = indexArray[2];
            Index d = indexArray[3];
            Coord center = (x[a]+x[b]+x[c]+x[d])*0.125;
            Coord pa = (x[a]+center)*(Real)0.666667;
            Coord pb = (x[b]+center)*(Real)0.666667;
            Coord pc = (x[c]+center)*(Real)0.666667;
            Coord pd = (x[d]+center)*(Real)0.666667;

            points[0].push_back(pa);
            points[0].push_back(pb);
            points[0].push_back(pc);

            points[1].push_back(pb);
            points[1].push_back(pc);
            points[1].push_back(pd);

            points[2].push_back(pc);
            points[2].push_back(pd);
            points[2].push_back(pa);

            points[3].push_back(pd);
            points[3].push_back(pa);
            points[3].push_back(pb);

            if(d_drawDirection.getValue() and d_anisotropyDirection.getValue().size() == _topology->getNbTetrahedra())
            {
                // msg_info() << "----- d_anisotropyDirection ----- " << i;
                Coord n = d_anisotropyDirection.getValue()[i];
                n/=n.norm();

                sofa::defaulttype::SolidTypes<double>::Transform t;
                t.getOrigin() = Vec3(0,0,0);
                Vec3 temp= n; // t.getOrientation().rotate(n);

                Coord tetraBarycenter = (x[a] + x[b] + x[c] + x[d]) / 4;
//                msg_info() << '[' << tetraBaryce  nter[0] << ',' << tetraBarycenter[1] << ',' << tetraBarycenter[2] << "],";
                Coord edge;
                Real tetraEdgeLenght;
                edge = x[a] - x[b];
                tetraEdgeLenght += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));
                edge = x[a] - x[c];
                tetraEdgeLenght += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));
                edge = x[a] - x[d];
                tetraEdgeLenght += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));
                edge = x[b] - x[c];
                tetraEdgeLenght += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));
                edge = x[b] - x[d];
                tetraEdgeLenght += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));
                edge = x[c] - x[c];
                tetraEdgeLenght += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));

                tetraEdgeLenght /= 6;
                Real tetraInsphere =  tetraEdgeLenght / sqrt(24);
                temp *= tetraInsphere ;
                type::Quat<SReal> q;
                q.fromMatrix(tetrahedronInf[i].rotation);
                q.normalize();
                temp = q.inverseRotate(temp);

                glLineWidth(3);
                vparams->drawTool()->drawLine(tetraBarycenter-temp,tetraBarycenter+temp, type::RGBAColor(0.,0.,0.,.8));
            }
            if(d_drawHeterogeneousTetra.getValue()) {
                // msg_info() << "----- d_drawHeterogeneousTetra ----- " << i;
                young_1 = d_youngModulus.getValue()[i];
                young_2 = d_anisotropyParameter.getValue()[i][1];
                if(young_1 - young_2 != 0)
                {
                    float col = (float)(abs(young_1 - young_2)/(young_1 + young_2));
                    //msg_warning() << col;
                    float fac = col * 0.5f;
                    type::RGBAColor color1(col      , 0.0f - fac , 1.0f-col,transparency);
                    type::RGBAColor color2(col      , 0.5f - fac , 1.0f-col,transparency);
                    type::RGBAColor color3(col      , 1.0f - fac , 1.0f-col,transparency);
                    type::RGBAColor color4(col+0.5f , 1.0f - fac , 1.0f-col,transparency);

                    vparams->drawTool()->drawTriangles(points[0],color2 );
                    vparams->drawTool()->drawTriangles(points[1],color2 );
                    vparams->drawTool()->drawTriangles(points[2],color2 );
                    vparams->drawTool()->drawTriangles(points[3],color2 );
                }
            }
            //            else if (anisotropyType == 4)
            //            {
            //                //msg_warning() << d_youngModulus.getValue()[i] << "  " << d_anisotropyParameter.getValue()[i][1];

            //                vparams->drawTool()->drawTriangles(points[0], defaulttype::Vec<4,float>(0.0,0.0,1.0,1.0));
            //                vparams->drawTool()->drawTriangles(points[1], defaulttype::Vec<4,float>(0.0,0.5,1.0,1.0));
            //                vparams->drawTool()->drawTriangles(points[2], defaulttype::Vec<4,float>(0.0,1.0,1.0,1.0));
            //                vparams->drawTool()->drawTriangles(points[3], defaulttype::Vec<4,float>(0.5,1.0,1.0,1.0));
            //            }
            for(unsigned int i=0 ; i<4 ; i++) points[i].clear();

            bool drawControlPoints = true;
            if(drawControlPoints)
            {
                const type::vector<Vec5>&  controlPoints = d_controlPoints.getValue();
                if (controlPoints.size() > 0)
                {
                    Coord CP;
                    for (size_t j=0;j<controlPoints.size();j++)
                    {
                        vparams->drawTool()->drawSphere(x[controlPoints[j][0]], 1.5,type::RGBAColor(102./255.,0,102./255.,1.0));//type::RGBAColor(230./255.,32./255.,32./255.,.8));
                    }
                }

            }

        }
        tetrahedronInfo.endEdit();
    }
#endif /* SOFA_NO_OPENGL */

}


}
