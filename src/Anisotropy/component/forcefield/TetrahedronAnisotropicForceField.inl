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
}


template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);

    this->Inherited::init();


    if (d_anisotropy.getValue() == "transverseIsotropic")
        elasticitySymmetry = TRANSVERSE_ISOTROPIC;
    else if (d_anisotropy.getValue() == "cubic")
        elasticitySymmetry = CUBIC;
    else {
        if (d_anisotropy.getValue() != "isotropic")
            msg_warning() << "Unknown elasticitySymmetry '" << d_anisotropy.getValue() << "', falling back to isotropic.";
        elasticitySymmetry = ISOTROPIC;
    }

    /// Init of tetrahedronInf which is a container which has, among other, mechanical info for each tetra:
    /// shapeVector, rotation, stiffnessVector ...
    _topology = this->getContext()->getMeshTopology();
    const std::vector< Tetrahedron > &tetrahedronArray = this->_topology->getTetrahedra();
    type::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    tetrahedronInf.resize(_topology->getNbTetrahedra());

    /// We are putting an elasticity tensor (E) for each tetra of our topology.
    /// It is decomposed via Kelvin modes into eigenvalues + eigentensors:
    ///   E = sum([k=1->N_modes] of lambda_k * P_k)   (Germain 2013, Eq 4.3)
    /// Two valid configurations:
    ///   - controlPoints set: IDW interpolates per-tetra params (direction, E_T, E_L)
    ///   - anisotropyDirections set with one entry per tetra: uniform or pre-assigned params
    const bool hasControlPoints = d_controlPoints.isSet() && d_controlPoints.getValue().size() > 0;
    const bool hasDirections    = d_anisotropyDirection.getValue().size() == _topology->getNbTetrahedra();

    if (hasControlPoints || hasDirections)
    {
        const size_t nbTetras = _topology->getNbTetrahedra();
        for (size_t i=0; i<nbTetras; ++i)
        {
            if (hasControlPoints)
                setMechanicalParametersFromControlPoints(i, tetrahedronArray[i]);

            if (i == 0)
            {
                if (d_anisotropyDirection.getValue().size() != nbTetras ||
                    d_anisotropyParameter.getValue().size() != nbTetras ||
                    d_youngModulus.getValue().size()         != nbTetras ||
                    d_poissonRatio.getValue().size()         != nbTetras)
                {
                    msg_error() << "Size mismatch before computeKelvinModesForElts:"
                                << " anisotropyDirection=" << d_anisotropyDirection.getValue().size()
                                << ", anisotropyParameter=" << d_anisotropyParameter.getValue().size()
                                << ", youngModulus=" << d_youngModulus.getValue().size()
                                << ", poissonRatio=" << d_poissonRatio.getValue().size()
                                << ", expected=" << nbTetras;
                    tetrahedronInfo.endEdit();
                    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
                    return;
                }
            }

            computeKelvinModesForElts(i, tetrahedronInf[i].eigenTensors, tetrahedronInf[i].eigenValues);
            initStiffnessVector(tetrahedronInf[i], tetrahedronArray[i]);
        }
        verifyStability();
    }
    else if (elasticitySymmetry == ISOTROPIC)
    {
        // No fiber direction needed — stiffness is built purely from Lamé coefficients (lambda, mu).
        updateLameCoefficients();
        for (size_t i = 0; i < _topology->getNbTetrahedra(); ++i)
            initStiffnessVector(tetrahedronInf[i], tetrahedronArray[i]);
    }
    else
        msg_error() << "Neither anisotropyDirections (size=" << d_anisotropyDirection.getValue().size()
                    << ", expected " << _topology->getNbTetrahedra() << ") nor controlPoints are set.";


    tetrahedronInfo.endEdit();
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    updateTopologyInfo=true;
}


template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::reinit()
{
    const std::vector<Tetrahedron>& tetrahedronArray = this->_topology->getTetrahedra();
    type::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());

    for (size_t i = 0; i < _topology->getNbTetrahedra(); ++i)
    {
        if (d_controlPoints.getValue().size() > 0)
            setMechanicalParametersFromControlPoints(i, tetrahedronArray[i]);
        computeKelvinModesForElts(i, tetrahedronInf[i].eigenTensors, tetrahedronInf[i].eigenValues);
        updateStiffnessVector(tetrahedronArray[i], tetrahedronInf[i]);
    }

    tetrahedronInfo.endEdit();
    verifyStability();
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

/// IDW interpolation of E_T, E_L, and fiber direction for one tetrahedron from d_controlPoints.
/// Weights = 1/dist^3 normalized over the d_IDWDepth nearest CPs (see generateListOfNorm).
/// ν_TL and G_L are NOT interpolated — Vec5 control points have no room for them (see seed block below).
template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::setMechanicalParametersFromControlPoints(size_t eltIndex,Tetra indexArray)
{
    AnisotropyDirectionArray &dir = *d_anisotropyDirection.beginEdit();
    type::vector<ParameterArray> &param = *d_anisotropyParameter.beginEdit();
    type::vector<Real> &youngModulus = *d_youngModulus.beginEdit();
    type::vector<Real> &poissonRatio = *d_poissonRatio.beginEdit();

    if (updateTopologyInfo and eltIndex == 0)
    {
        // Save user-provided scalar values before vectors are resized.
        // IDW only interpolates E_T, E_L, and direction (Vec5 has no room for ν_TL or G_L).
        // ν_T  comes from d_poissonRatio[0] (user-supplied, e.g. via Python poissonRatio=…).
        // ν_TL and G_L come from d_anisotropyParameters[0][2..3] (user-supplied, or defaults below).
        // To support spatially varying ν_TL/G_L, extend controlPoints to Vec7.
        const Real defaultNuT  = poissonRatio.size() >= 1 ? poissonRatio[0] : Real(0.4161);
        const ParameterArray defaultAnisoPar = param.size() >= 1
                                               ? param[0]
                                               : ParameterArray{2, 0., Real(0.4161), Real(0.)};
        dir.resize(_topology->getNbTetrahedra());
        param.resize(_topology->getNbTetrahedra());
        youngModulus.resize(_topology->getNbTetrahedra());
        poissonRatio.resize(_topology->getNbTetrahedra());
        for (size_t i=0; i<_topology->getNbTetrahedra(); ++i)
        {
            param[i] = {defaultAnisoPar[0], 0., defaultAnisoPar[2], defaultAnisoPar[3]};
            poissonRatio[i] = defaultNuT;
        }
    }

    const VecCoord& x = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();
    const Coord tetraBarycenter= (x[indexArray[0]] + x[indexArray[1]] + x[indexArray[2]] + x[indexArray[3]]) / 4;

    type::vector<std::pair<Real, int> > listOfNorm = generateListOfNorm(tetraBarycenter);

    // controlPoints[j] = [nodeIndex, E_T, E_L, θ, φ]
    IDWdata(listOfNorm, youngModulus[eltIndex], 1); // E_T (Young Modulus Transverse)
    IDWdata(listOfNorm, param[eltIndex][1],     2); // E_L (Young Modulus Longitudinal)
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

    const VecCoord& x = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();

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

/// Computes per-element Kelvin eigenvalues (λ_k) and eigentensors (N_k) of the elasticity tensor.
/// The elasticity tensor decomposes as E = Σ_k λ_k (N_k ⊗ N_k) — Germain 2013, Eq. 4.3.
/// CUBIC:                3 eigenvalues, 6 tensors (Nd, Ne, Np, Ns1-3)         — Germain §4.4.1, Eqs 4.39-4.41
/// TRANSVERSE_ISOTROPIC: 4 eigenvalues, 6 tensors (Nh1, Np, Ns3, Nh2, Ns1, Ns2) — Germain §4.4.7, Eqs 4.118-4.122
/// The local frame (n = fiber axis, v1/v2 = transverse orthonormal basis) is built from d_anisotropyDirection[eltIndex].
template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::computeKelvinModesForElts(size_t eltIndex,EigenTensors &eigenTensors , EigenValues &eigenValues)
{
    if (elasticitySymmetry != ISOTROPIC)
    {
        Vec4 anisotropyParameter=d_anisotropyParameter.getValue()[eltIndex];
        // Vec10 anisotropyParameter=d_anisotropyParameter.getValue()[eltIndex];
        Real youngModulus = d_youngModulus.getValue()[eltIndex];
        Real poissonRatio = d_poissonRatio.getValue()[eltIndex];

        // n: fiber/symmetry axis (e3 in Germain's notation); v1, v2: orthonormal transverse basis
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

            // Voigt stiffness entries from isotropic E and ν, modified by the Zener anisotropy ratio A.
            // A=1 → isotropic (c44 = μ), A≠1 → stiffer (A>1) or softer (A<1) along cube-face diagonals.
            Real c11 = youngModulus * (1 - poissonRatio) / (1 - poissonRatio - 2 * poissonRatio * poissonRatio);
            Real c12 = youngModulus * poissonRatio / (1 - poissonRatio - 2 * poissonRatio * poissonRatio);
            Real c44 = anisotropyRatio * (c11 - c12) / 2;

            // Germain §4.4.1, Eqs 4.39-4.41: 3 eigenvalues with multiplicities 1, 2, 3
            Real eigen1 = c11 + 2 * c12;    // (4.39) bulk — dilatation mode Nd
            Real eigen2 = c11 - c12;        // (4.40) isotropic shear — modes Ne, Np
            Real eigen3 = 2 * c44;          // (4.41) cubic shear — modes Ns1, Ns2, Ns3

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


            eigenTensors[0] = Nd;  eigenValues[0] = eigen1;
            eigenTensors[1] = Ne;  eigenValues[1] = eigen2;
            eigenTensors[2] = Np;  eigenValues[2] = eigen2;
            eigenTensors[3] = Ns1; eigenValues[3] = eigen3;
            eigenTensors[4] = Ns2; eigenValues[4] = eigen3;
            eigenTensors[5] = Ns3; eigenValues[5] = eigen3;

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
            // get the constants from the young modulus, Poisson ratio and anisotropy ratio.
            Real youngModulusTransverse = youngModulus;
            Real poissonRatioTransverse = poissonRatio;
            Real youngModulusLongitudinal=anisotropyParameter[1];
            Real poissonRatioTransverseLongitudinal=anisotropyParameter[2];
            Real shearModulusLongitudinal=anisotropyParameter[3];

            Real poissonRatioLongitudinalTransverse= poissonRatioTransverseLongitudinal*youngModulusLongitudinal/youngModulusTransverse;

            // Voigt stiffness entries (fiber axis = e3): c11=c22 (transverse), c33 (longitudinal),
            // c12 (ν_T coupling), c13=c23 (mixed), c44=c55=G_L (fiber shear), c55=(c11-c12)/2=G_T (transverse shear)
            Real gammaInv = 1 - pow(poissonRatioTransverse,2)
                              - 2*poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal
                              - 2*poissonRatioTransverse*poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal;
            // Stability is checked per-element in the calling loop (init/reinit) to aggregate counts.
            Real gamma=1/gammaInv;

            Real c11=youngModulusTransverse*gamma*(1-poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal); // = c22
            Real c33=youngModulusLongitudinal*gamma*(1-poissonRatioTransverse*poissonRatioTransverse);
            Real c44= youngModulusTransverse/(2*(1+poissonRatioTransverse)); // G_T (transverse shear) — see eigen4 bug below
            Real c55= c44; // = G_T (same as c44 for this implementation)
            Real c12=youngModulusTransverse*gamma*(poissonRatioTransverse+poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);
            Real c13=youngModulusTransverse*gamma*(poissonRatioLongitudinalTransverse+poissonRatioTransverse*poissonRatioTransverseLongitudinal);


            ///////////////////////////////////////////////////////////////////////////////////////////////////



            // (4.119)
            Real talpha=sqrt(2.0f)*(c11+c12-c33)/(4*c13);
            Real alpha=atan(talpha);
            Real salpha=sin(alpha);
            Real calpha=cos(alpha);
            Real secalpha=1/calpha;

            // Germain §4.4.7, Eqs 4.118-4.119 (hexagonal tensor, e3 axis)
            Real eigen1 = c33 + M_SQRT2 * c13 * (talpha + secalpha); // λ1: mixed bulk/fiber dilatation
            Real eigen2 = c11 - c12;                                  // λ2: transverse shear 2G_T (modes Np, Ns3)
            Real eigen3 = c33 + M_SQRT2 * c13 * (talpha - secalpha); // λ3: mixed extension/compression
            Real eigen4 = 2 * shearModulusLongitudinal; // 2  * c44; // BUG: should be 2*shearModulusLongitudinal (2G_L). See project_known_bugs.md.
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
                if (eltIndex == 0) 
                    msg_info() << "---------> \nTRANSVERSE_ISOTROPIC but set with ISOTROPIC mechanical parameter switch to ISOTROPIC kelvin mode decomposition/projection";
                v1 = Coord(1,0,0);
                v2 = Coord(0,1,0);
                n  = Coord(0,0,1);

                // (4.29)
                eigen1 = c11 + 2*c12;
                eigen2 = c11 - c12;
                eigen3 = eigen2;
                eigen4 = eigen2;

                Nh1.identity(); // == Nd
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
}

template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::verifyStability() const
{
    const type::vector<TetrahedronRestInformation>& tetrahedronInf = tetrahedronInfo.getValue();
    const size_t nbTetras = tetrahedronInf.size();
    if (nbTetras == 0 || elasticitySymmetry == ISOTROPIC) return;

    size_t nbUnstable = 0, firstUnstable = size_t(-1);

    for (size_t i = 0; i < nbTetras; ++i)
    {
        const EigenValues& ev = tetrahedronInf[i].eigenValues;
        Real minEig = std::numeric_limits<Real>::max();
        Real maxEig = std::numeric_limits<Real>::lowest();
        bool eltUnstable = false;
        for (size_t k = 0; k < ev.size(); ++k) {
            if (ev[k] < minEig) minEig = ev[k];
            if (ev[k] > maxEig) maxEig = ev[k];
            if (ev[k] <= Real(0)) eltUnstable = true;
        }

        if (eltUnstable) {
            if (nbUnstable == 0) {
                // Detailed diagnosis for the first unstable element only.
                const Vec4 ap  = d_anisotropyParameter.getValue()[i];
                Real E_T       = d_youngModulus.getValue()[i];
                Real nu_T      = d_poissonRatio.getValue()[i];
                Real E_L       = ap[1];
                Real nu_TL     = ap[2];
                Real nu_LT     = nu_TL * E_L / E_T;
                Real gammaInv  = 1 - nu_T*nu_T - 2*nu_LT*nu_TL - 2*nu_T*nu_LT*nu_TL;
                Real nuTLBound = (E_L > Real(0)) ? std::sqrt((1-nu_T) / (2*E_L/E_T)) : Real(0);
                msg_error() << "TI UNSTABLE at elt " << i << " (first occurrence):"
                            << "\n  min Kelvin eigenvalue = " << minEig
                            << "\n  gammaInv = " << gammaInv << (gammaInv <= 0 ? " <= 0 (stiffness tensor is indefinite)" : "")
                            << "\n  E_L/E_T = " << E_L/E_T
                            << "\n  nu_TL = " << nu_TL << "  stability bound: |nu_TL| < " << nuTLBound
                            << "\n  nu_LT (derived) = " << nu_LT
                            << "\n  Use tools/ti_material_stability.py to find valid parameter ranges.";
            }
            ++nbUnstable;
            if (firstUnstable == size_t(-1)) firstUnstable = i;
        } else if (nbUnstable == 0 && i == 0 && minEig > Real(0)) {
            // Condition number check on first element if stable.
            Real cond = maxEig / minEig;
            if (cond > Real(1e4))
                msg_warning() << "TI stiffness condition number ~" << cond
                              << " at elt 0 (max/min Kelvin eigenvalue)."
                              << " Iterative solvers (CGLinearSolver) may converge slowly.";
        }
    }

    if (nbUnstable > 0)
        msg_error() << nbUnstable << "/" << nbTetras
                    << " elements have at least one negative Kelvin eigenvalue."
                    << " First unstable element: " << firstUnstable
                    << ". Simulation WILL diverge — fix material parameters.";
}

template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::initStiffnessVector(TetrahedronRestInformation &my_tinfo,
                                                                      const Tetrahedron &tt)
{
    // set array to zero
    std::fill(my_tinfo.stiffnessVector.begin(), my_tinfo.stiffnessVector.end(), Mat3x3());


    typename DataTypes::Coord point[4];
    const typename DataTypes::VecCoord &restPosition = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();

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
}

template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::updateStiffnessVector(const Tetrahedron &tetra,TetrahedronRestInformation &my_tinfo)
{
    const typename DataTypes::VecCoord &position = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();

    // set array to zero
    std::fill(my_tinfo.stiffnessVector.begin(), my_tinfo.stiffnessVector.end(), Mat3x3());

    typename DataTypes::Coord point[4];
    // store the point position
    for (size_t j = 0; j < 4; ++j)
        point[j] = (position)[tetra[j]];

    computeTetrahedronStiffnessEdgeMatrixForElts(my_tinfo.stiffnessVector,my_tinfo.eigenTensors,my_tinfo.eigenValues,point);
}

/// Assembles the 6 edge stiffness 3×3 matrices K[0..5] for one tetrahedron.
/// Shape vectors b_j = ±cross(opposite edges)/vol are the P1 barycentric coordinate gradients (const per tet).
/// ISOTROPIC:    K[j]_mn = λ b_k[n]b_l[m] + μ(b_l[n]b_k[m] + (b_l·b_k)δ_mn)   — Delingette 2019, Prop 1
/// ANISOTROPIC:  K[j] = (|vol|/6) Σ_i λ_i A_i (b_l ⊗ b_k) A_i                 — Delingette 2019, Prop 3
/// Edge ordering j: {(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)} — see edgesInTetrahedronArray.
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

            // K[j] += V_K * Σ_i λ_i A_i (b_l⊗b_k) A_i — Delingette 2019, Proposition 3
            // tmp = b_l⊗b_k (dyadic product of shape vectors for edge j=(k,l))
            Mat3x3 tmp=type::dyad(shapeVector[l],shapeVector[k]);
            for(i=0;i<6;++i) {
                stiffnessVector[j]+=(eigenValues[i]*eigenTensors[i]*tmp*eigenTensors[i])*fabs(volume)/6;
            }
        }
    }
}


/// Extracts the rigid-body rotation from tetrahedron edge vectors dp[0..5] via QR decomposition.
/// Builds an orthonormal frame stored row-wise: row0=e_x (along dp[0]), row2=normalize(e_x × dp[1]), row1=row2 × row0.
/// Used at rest to compute restRotation and each frame to compute the current frame S; R = S^T * restRotation.
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

/// Corotational linear elasticity: per-tetra rigid rotation R is extracted via QR decomposition each step.
/// Forces are computed in the element rest frame then rotated back to world frame: f = R^T * (K * R * d).
/// Also caches rotatedStiffnessVector = R * K * R^T per edge for use by addDForce (implicit solvers).
template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */,
                                                                        DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & /*dataV*/ )
{
    SOFA_UNUSED(mparams);

    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& x  = dataX.getValue();
    const VecCoord& x0= this->mstate->read(core::vec_id::read_access::restPosition)->getValue();
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

    // Debug accumulators (reset each addForce call)
    Real dbg_maxRestForce = 0, dbg_maxDpos = 0, dbg_maxRerr = 0;
    size_t dbg_worstForceTet = 0, dbg_worstDposTet = 0, dbg_nanTet = size_t(-1);

    Coord dpos,sv;
    typename VecElement::const_iterator it;
    for(it=_indexedElements->begin(), i = 0 ; it!=_indexedElements->end(); ++it,++i)
    {
        Tetra index = *it;
        tetinfo=&tetrahedronInf[i];
        size_t nbPoints=index.size();
        const  TetraEdgesStiffness &stiffnessArray=getStiffnessArray(tetinfo);


        Mat3x3 deformationGradient,S,R;
        Coord dpp[6];
        for (j=0; j<6; ++j)
        {
            dpp[j]=x[tetinfo->v[edgesInTetrahedronArray[j][1]]]-x[tetinfo->v[edgesInTetrahedronArray[j][0]]];
        }

        // S: QR frame of current edges; R = S^T * restRotation maps world-frame vectors to rest-frame
        computeQRRotation(S,dpp);
        R=S.transposed()*tetinfo->restRotation;

        // Debug: R orthogonality residual (cheap: 9 multiplies)
        if (this->f_printLog.getValue()) {
            Mat3x3 RtR = R.transposed() * R;
            Real Rerr = 0;
            for (int a=0;a<3;a++) for (int b=0;b<3;b++) {
                Real d = RtR[a][b] - (a==b ? Real(1) : Real(0));
                Rerr = std::max(Rerr, std::abs(d));
            }
            if (Rerr > dbg_maxRerr) dbg_maxRerr = Rerr;
        }

        // tetinfo->rotation = R^T: maps rest-frame vectors back to world frame (used for force back-rotation)
        tetinfo->rotation=R.transposed();
        std::fill(force.begin(),force.end(),Coord());

        tetinfo->rotatedStiffnessVector.clear();

        // loop over each entry in the stiffness vector
        for (l=0,j=0; j<nbPoints; ++j) {
            v0 = index[j];
            for ( k=j+1; k<nbPoints; ++k,++l) {
                v1 = index[k];
                dpos=x[v0]-x[v1];
                // d_rest = R^T * (x_v0 - x_v1) - (x0_v0 - x0_v1): displacement mapped to rest frame
                dpos=tetinfo->rotation*dpos-(x0[v0]-x0[v1]);

                // Debug: track max strain magnitude
                {
                    Real dn = dpos.norm();
                    if (!std::isfinite(dn)) {
                        if (dbg_nanTet == size_t(-1)) dbg_nanTet = i;
                    } else if (dn > dbg_maxDpos) { dbg_maxDpos = dn; dbg_worstDposTet = i; }
                }

                // forces accumulated in rest frame (Newton's 3rd law for the two endpoints)
                force[k]-=stiffnessArray[l]*dpos;
                force[j]+=stiffnessArray[l].multTranspose(dpos);

                // cache K_rot = R * K * R^T (world-frame stiffness) — unused by addDForce/addKToMatrix, kept for debugging
                Mat3x3 mat=R*stiffnessArray[l]*tetinfo->rotation;
                tetinfo->rotatedStiffnessVector[l] = mat;

            }
        }

        // Debug: track max rest-frame force magnitude
        for (size_t n=0; n<nbPoints; ++n) {
            Real fn = force[n].norm();
            if (!std::isfinite(fn)) {
                if (dbg_nanTet == size_t(-1)) dbg_nanTet = i;
            } else if (fn > dbg_maxRestForce) { dbg_maxRestForce = fn; dbg_worstForceTet = i; }
        }

        for (j=0; j<nbPoints; ++j)
        {
            f[index[j]]+=R*force[j];
        }

    }

    // Debug summary: print every call when printLog=true, or on any NaN, or when forces are suspiciously large
    if (this->f_printLog.getValue() || dbg_nanTet != size_t(-1) || dbg_maxRestForce > Real(1e6)) {
        msg_info() << "[addForce] maxRestForce=" << dbg_maxRestForce << " @tet=" << dbg_worstForceTet
                   << "  maxDpos(strain)=" << dbg_maxDpos << " @tet=" << dbg_worstDposTet
                   << "  maxRerr(R orth)=" << dbg_maxRerr;
        if (dbg_nanTet != size_t(-1))
            msg_warning() << "[addForce] NaN/Inf detected at tet=" << dbg_nanTet
                          << " — positions or stiffness may be corrupt";
    }

    updateMatrix=true; // useless normally
    tetrahedronInfo.endEdit();
    dataF.endEdit();
}


/// Computes df = kFactor * K_rot * dx for implicit time integration.
/// Recomputes R from current positions each call so K_rot = R*K*R^T is always fresh.
/// This is necessary for TI materials: unlike isotropic K, K_rot is NOT rotation-invariant,
/// so using a cached R from a previous addForce call produces wrong tangent stiffness under
/// fast rigid-body motion, causing the solver to diverge.
/// kFactor = (h² * stiffness_scale + h * rayleighStiffness) incorporates Rayleigh damping.
template <class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX )
{
    VecDeriv& df       = *(datadF.beginEdit());
    const VecCoord& dx =   datadX.getValue()  ;
    Real kFactor = sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams,this->rayleighStiffness.getValue());
    const VecCoord& x  = this->mstate->read(core::vec_id::read_access::position)->getValue();
    type::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());

    Coord dpos;
    size_t i,j,k,v0,v1,rank;

    TetrahedronRestInformation *tetinfo;
    const VecElement * _indexedElements = & (_topology->getTetrahedra());

    Real dbg_maxDf = 0, dbg_maxDx = 0;
    size_t dbg_worstDfTet = 0, dbg_worstDxTet = 0, dbg_nanTet = size_t(-1);

    typename VecElement::const_iterator it;
    for(it=_indexedElements->begin(), i = 0 ; it!=_indexedElements->end(); ++it,++i)
    {
        Tetra index = *it;
        const size_t nbPoints = index.size();
        tetinfo=&tetrahedronInf[i];

        // Recompute R from current positions (not cached from last addForce).
        Coord dpp[6];
        for (j=0; j<6; ++j)
            dpp[j] = x[tetinfo->v[edgesInTetrahedronArray[j][1]]] - x[tetinfo->v[edgesInTetrahedronArray[j][0]]];
        Mat3x3 S, R;
        computeQRRotation(S, dpp);
        R = S.transposed() * tetinfo->restRotation;
        const Mat3x3 Q = R.transposed(); // Q = R^T: world → rest frame

        sofa::type::vector<Deriv> dforce;
        dforce.resize(nbPoints);

        for (rank=0,j=0; j<nbPoints; ++j)
        {
            v0 = index[j];
            for ( k=j+1; k<nbPoints; ++k,++rank)
            {
                v1 = index[k];
                dpos = Q * (dx[v0]-dx[v1]);          // rotate dx to rest frame

                // Debug: track max |dx| input
                {
                    Real dn = dpos.norm();
                    if (!std::isfinite(dn)) { if (dbg_nanTet == size_t(-1)) dbg_nanTet = i; }
                    else if (dn > dbg_maxDx) { dbg_maxDx = dn; dbg_worstDxTet = i; }
                }

                const Coord fk = tetinfo->stiffnessVector[rank] * dpos;
                const Coord fj = tetinfo->stiffnessVector[rank].multTranspose(dpos);
                dforce[k] -= R * fk * kFactor;       // rotate force back to world frame
                dforce[j] += R * fj * kFactor;
            }
        }

        // Debug: track max |df| output
        for (size_t n=0; n<nbPoints; ++n) {
            Real fn = dforce[n].norm();
            if (!std::isfinite(fn)) { if (dbg_nanTet == size_t(-1)) dbg_nanTet = i; }
            else if (fn > dbg_maxDf) { dbg_maxDf = fn; dbg_worstDfTet = i; }
        }

        for (j=0; j<nbPoints; ++j)
        {
            df[index[j]]+=dforce[j];
        }
    }

    if (this->f_printLog.getValue() || dbg_nanTet != size_t(-1) || dbg_maxDf > Real(1e6)) {
        msg_info() << "[addDForce] kFactor=" << kFactor
                   << "  maxDf=" << dbg_maxDf << " @tet=" << dbg_worstDfTet
                   << "  maxDx(rest)=" << dbg_maxDx << " @tet=" << dbg_worstDxTet;
        if (dbg_nanTet != size_t(-1))
            msg_warning() << "[addDForce] NaN/Inf detected at tet=" << dbg_nanTet;
    }

    datadF.endEdit();
    tetrahedronInfo.endEdit();
}


/// Assembles the global tangent stiffness matrix for direct solvers (e.g. SparseLDLSolver).
/// Recomputes R from current positions so K_rot = R*K*R^T matches the actual configuration.
/// Adds 4 blocks per edge pair (j,l): K at (l,l), -K at (l,j), -K^T at (j,l), K^T at (j,j).
template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::addKToMatrix(sofa::linearalgebra::BaseMatrix *mat, SReal k, unsigned int &offset)
{
    const VecCoord& x  = this->mstate->read(core::vec_id::read_access::position)->getValue();
    const type::vector<TetrahedronRestInformation> tetrahedronInf = *(tetrahedronInfo.beginEdit());
    size_t i,j,l,rank,n,m;

    const TetrahedronRestInformation *tetinfo;
    const VecElement * _indexedElements = & (_topology->getTetrahedra());

    Real dbg_maxKrot = 0;
    size_t dbg_worstKrotTet = 0, dbg_nanTet = size_t(-1);

    typename VecElement::const_iterator it;
    for(it=_indexedElements->begin(), i = 0 ; it!=_indexedElements->end(); ++it,++i)
    {
        Tetra index = *it;
        const size_t nbPoints = index.size();
        tetinfo=&tetrahedronInf[i];

        // Recompute R from current positions.
        Coord dpp[6];
        for (j=0; j<6; ++j)
            dpp[j] = x[tetinfo->v[edgesInTetrahedronArray[j][1]]] - x[tetinfo->v[edgesInTetrahedronArray[j][0]]];
        Mat3x3 S, R;
        computeQRRotation(S, dpp);
        R = S.transposed() * tetinfo->restRotation;
        const Mat3x3 Q = R.transposed();

        // Build fresh K_rot[rank] = R * K[rank] * R^T for all 6 edge pairs.
        TetraEdgesStiffness rotatedK;
        for (rank=0; rank<6; ++rank) {
            rotatedK[rank] = R * tetinfo->stiffnessVector[rank] * Q;

            // Debug: track max absolute entry in K_rot
            for (n=0; n<3; ++n) for (m=0; m<3; ++m) {
                Real v = std::abs(rotatedK[rank][n][m]);
                if (!std::isfinite(v)) { if (dbg_nanTet == size_t(-1)) dbg_nanTet = i; }
                else if (v > dbg_maxKrot) { dbg_maxKrot = v; dbg_worstKrotTet = i; }
            }
        }

        for (rank=0,j=0; j<nbPoints; ++j) {
            for ( l=j+1; l<nbPoints; ++l,++rank) {
                for (n=0; n<3; ++n) {
                    for (m=0; m<3; ++m) {
                        Mat3x3 rotatedKTranspose = rotatedK[rank];
                        rotatedKTranspose.transpose();
                        mat->add(3*index[l]+ n + offset, 3*index[l]+ m + offset,   rotatedK[rank][n][m] * k);
                        mat->add(3*index[l]+ n + offset, 3*index[j]+ m + offset, - rotatedK[rank][n][m] * k);
                        mat->add(3*index[j]+ n + offset, 3*index[l]+ m + offset, - rotatedKTranspose[n][m] * k);
                        mat->add(3*index[j]+ n + offset, 3*index[j]+ m + offset,   rotatedKTranspose[n][m] * k);
                    }
                }
            }
        }
    }

    if (this->f_printLog.getValue() || dbg_nanTet != size_t(-1) || dbg_maxKrot > Real(1e9)) {
        msg_info() << "[addKToMatrix] k=" << k
                   << "  maxKrot=" << dbg_maxKrot << " @tet=" << dbg_worstKrotTet;
        if (dbg_nanTet != size_t(-1))
            msg_warning() << "[addKToMatrix] NaN/Inf in K_rot at tet=" << dbg_nanTet;
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
SReal TetrahedronAnisotropicForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord& dataX) const
{
    const VecCoord& x  = dataX.getValue();
    const VecCoord& x0 = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();
    const type::vector<TetrahedronRestInformation>& tetrahedronInf = tetrahedronInfo.getValue();
    const VecElement* elements = &(_topology->getTetrahedra());

    SReal energy = 0;
    size_t i = 0;
    for (auto it = elements->begin(); it != elements->end(); ++it, ++i)
    {
        const Tetra& index = *it;
        const TetrahedronRestInformation& tetinfo = tetrahedronInf[i];

        // Rest-frame nodal displacements: u_j = R^T * x_j - x0_j
        // tetinfo.rotation stores R^T (world->rest), cached by the last addForce call.
        Coord u[4];
        for (size_t j = 0; j < 4; ++j)
            u[j] = tetinfo.rotation * x[index[j]] - x0[index[j]];

        // W = ½ Σ_{edges (j,k)} [ (K^T·d)·u_j − (K·d)·u_k ]
        // where d = u_j - u_k.  Mirrors addForce: force[j] += K^T·d, force[k] -= K·d.
        // For symmetric K this reduces to ½ d^T·K·d; the full form is exact for any K.
        size_t rank = 0;
        for (size_t j = 0; j < 4; ++j)
            for (size_t k = j + 1; k < 4; ++k, ++rank) {
                const Coord d = u[j] - u[k];
                energy += SReal(0.5) * (dot(tetinfo.stiffnessVector[rank].multTranspose(d), u[j])
                                      - dot(tetinfo.stiffnessVector[rank] * d, u[k]));
            }
    }
    return energy;
}


template<class DataTypes>
void TetrahedronAnisotropicForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

#ifndef SOFA_NO_OPENGL
    if(this->d_componentState.getValue() == sofa::core::objectmodel::ComponentState::Invalid) return;
    if (!vparams->displayFlags().getShowForceFields()) return;
    {

        type::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
        const VecElement * _indexedElements = & (_topology->getTetrahedra());

        const VecCoord& x = this->mstate->read(core::vec_id::read_access::position)->getValue();

        size_t i;
        size_t nbTetrahedra=_topology->getNbTetrahedra();

        const bool hasAnisotropyParam = d_anisotropyParameter.getValue().size() == nbTetrahedra;
        Real young_1, young_2, anisotropyType;
        float transparency = d_transparency.getValue();
        typename VecElement::const_iterator it;

        for(it=_indexedElements->begin(), i = 0 ; it!=_indexedElements->end(); ++it,++i)
        {
            anisotropyType = hasAnisotropyParam ? d_anisotropyParameter.getValue()[i][0] : Real(0);


            Tetra indexArray = *it;
            std::vector< type::Vec3 > points[4];

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
                Coord n = d_anisotropyDirection.getValue()[i];
                n/=n.norm();

                sofa::defaulttype::SolidTypes<double>::Transform t;
                t.getOrigin() = Vec3(0,0,0);
                Vec3 temp= n; // t.getOrientation().rotate(n);

                Coord tetraBarycenter = (x[a] + x[b] + x[c] + x[d]) / 4;
                Coord edge;
                Real tetraEdgeLength = 0;
                edge = x[a] - x[b];
                tetraEdgeLength += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));
                edge = x[a] - x[c];
                tetraEdgeLength += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));
                edge = x[a] - x[d];
                tetraEdgeLength += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));
                edge = x[b] - x[c];
                tetraEdgeLength += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));
                edge = x[b] - x[d];
                tetraEdgeLength += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));
                edge = x[c] - x[d];
                tetraEdgeLength += sqrt(pow(edge[0], 2) + pow(edge[1], 2) + pow(edge[2], 2));

                tetraEdgeLength /= 6;
                Real tetraInsphere = tetraEdgeLength / sqrt(24);
                temp *= tetraInsphere ;
                type::Quat<SReal> q;
                q.fromMatrix(tetrahedronInf[i].rotation);
                q.normalize();
                temp = q.inverseRotate(temp);

                glLineWidth(3);
                vparams->drawTool()->drawLine(tetraBarycenter-temp,tetraBarycenter+temp, type::RGBAColor(0.,0.,0.,.8));
            }
            if(d_drawHeterogeneousTetra.getValue() && hasAnisotropyParam) {
                young_1 = d_youngModulus.getValue()[i];
                young_2 = d_anisotropyParameter.getValue()[i][1];
                if(young_1 - young_2 != 0)
                {
                    float col = (float)(abs(young_1 - young_2)/(young_1 + young_2));
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
            for(unsigned int i=0 ; i<4 ; i++) points[i].clear();

        }

        if (d_drawHeterogeneousTetra.getValue() && d_controlPoints.isSet() && d_controlPoints.getValue().size() > 0)
        {
            const type::vector<Vec5>& controlPoints = d_controlPoints.getValue();
            for (size_t j = 0; j < controlPoints.size(); ++j)
                vparams->drawTool()->drawSphere(x[controlPoints[j][0]], 1.5, type::RGBAColor(102./255., 0, 102./255., 1.0));
        }

        tetrahedronInfo.endEdit();
    }
#endif /* SOFA_NO_OPENGL */


}


}
